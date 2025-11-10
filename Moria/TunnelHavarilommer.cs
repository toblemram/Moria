using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using System.Drawing;
using System.IO;

namespace Tunnel.GH
{
    public class GH_TunnelHavarilommer : GH_Component
    {
        public GH_TunnelHavarilommer() : base(
            "Tunnel.Havarilommer",
            "Havarilommer",
            "Lager havarilommer etter N500 30–30–30-regel langs en tunnelbrep. " +
            "Profiltype styrer grunnbredde, nisjen gir +3,0 m ekstra bredde.",
            "Tunnel",
            "Details")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddBrepParameter(
                "Tunnel Brep",
                "B",
                "Tunnelbrep (fra Txx-profilkomponenten).",
                GH_ParamAccess.item);

            p.AddCurveParameter(
                "Path",
                "P",
                "Senterlinje/path som tunnelen er sveipet langs.",
                GH_ParamAccess.item);

            p.AddTextParameter(
                "Profiltype",
                "T",
                "Tunnelprofil, f.eks. T9,5 / T10,5 / T12,5 / T13 / T13,5 / T14.",
                GH_ParamAccess.item,
                "T14");

            p.AddNumberParameter(
                "Stasjoner høyre",
                "SR",
                "Midtstasjoner (m) langs path for havarilommer på høyre side (f.eks. 3, 60).",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "Stasjoner venstre",
                "SL",
                "Midtstasjoner (m) langs path for havarilommer på venstre side.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "Lommehøyde",
                "H",
                "Lommens høyde over bunn (m). Sett iht. N500.",
                GH_ParamAccess.item,
                2.4);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter(
                "Tunnel m/lommer",
                "B",
                "Tunnelbrep der havarilommer er slått sammen (BooleanUnion).",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Lommebreps",
                "Lommer",
                "Brep-er for hver genererte havarilomme (debug).",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "Lommesentre",
                "Ctr",
                "Senterpunkt (midt i volumet) for hver havarilomme.",
                GH_ParamAccess.list);

            p.AddTextParameter(
                "Info",
                "i",
                "Status og diagnostikk.",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "BadIntersectionPoints",
                "Bad",
                "Alle BadIntersectionPoints fra boolean-operasjoner (debug).",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "NakedEdgePoints",
                "Naked",
                "Alle NakedEdgePoints fra boolean-operasjoner (debug).",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "NonManifoldEdgePoints",
                "NonM",
                "Alle NonManifoldEdgePoints fra boolean-operasjoner (debug).",
                GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Brep tunnel = null;
            Curve path = null;
            string profileType = "T14";
            var stR = new List<double>();
            var stL = new List<double>();
            double lommeHoyde = 2.4;

            if (!da.GetData(0, ref tunnel) || tunnel == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler tunnelbrep.");
                return;
            }

            if (!da.GetData(1, ref path) || path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler path-kurve.");
                return;
            }

            da.GetData(2, ref profileType);
            da.GetDataList(3, stR);
            da.GetDataList(4, stL);
            da.GetData(5, ref lommeHoyde);

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;

            // Konstanter – juster ved behov
            const double NisjeExtra = 3.0;      // ekstra bredde i nisjen
            const double rampLen = 30.0;
            const double flatLen = 30.0;
            const double BottomLift = 0.05;     // løft bunn 5 cm for å unngå koplan snitt
            const double InnerInset = 0.20;     // lomme stikker 20 cm inn i tunnelen
            const int sampleCount = 17;
            const double tolTightFactor = 0.5;  // strammere tol ved retry
            const double tolLooseFactor = 2.0;  // slakkere tol ved retry

            var info = new List<string>();
            var lommer = new List<Brep>();
            var sentre = new List<Point3d>();
            var allBadPts = new List<Point3d>();
            var allNakedPts = new List<Point3d>();
            var allNonManifoldPts = new List<Point3d>();

            // ---- Profilbredde Bt → halv bredde ----
            profileType = profileType.Replace(",", ".").ToUpper();
            if (!ProfileWidths.TryGetValue(profileType, out double Bt))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Ugyldig profiltype '{profileType}'. Gyldige: {string.Join(", ", ProfileWidths.Keys)}");
                return;
            }

            double baseOffset = Bt / 2.0;    // skulder/vegg fra senterlinje

            double curveLength = path.GetLength();
            if (curveLength <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path-kurven er for kort/ugyldig.");
                return;
            }

            double totalLen = 2 * rampLen + flatLen; // 90 m
            double epsWidth = Math.Min(0.05, NisjeExtra * 0.02);

            // --- Hjelpefunksjon: gjør brep mest mulig solid, slik busslommer ofte gjør ---
            Brep PrepareSolid(Brep b, string label)
            {
                if (b == null) return null;

                if (!b.IsSolid)
                {
                    info.Add($"{label}: IsSolid == false – prøver CapPlanarHoles + Repair.");

                    var capped = b.CapPlanarHoles(tol);
                    if (capped != null)
                    {
                        info.Add($"{label}: CapPlanarHoles OK.");
                        b = capped;
                    }
                    else
                    {
                        info.Add($"{label}: CapPlanarHoles ga null – ingen planare hull?");
                    }

                    b.Repair(tol);
                    if (b.IsSolid)
                        info.Add($"{label}: Repair OK – brep er nå solid.");
                    else
                        info.Add($"{label}: Repair kjørt men brep er fortsatt ikke solid.");
                }
                else
                {
                    info.Add($"{label}: IsSolid == true.");
                }

                return b;
            }

            // ---------- Lokal funksjon: generer én nisje som Brep ----------
            Brep LagNisje(double sCenter, bool isRightSide)
            {
                double s0 = sCenter - (rampLen + flatLen / 2.0); // start
                double s1 = sCenter - (flatLen / 2.0);           // ramp ferdig
                double s2 = sCenter + (flatLen / 2.0);           // ramp starter ned
                double s3 = sCenter + (rampLen + flatLen / 2.0); // slutt

                if (s0 < 0 || s3 > curveLength)
                {
                    info.Add($"Nisje @ s={sCenter:0.###} m går utenfor path (0–{curveLength:0.###} m). Hopper over.");
                    return null;
                }

                var bi = new List<Point3d>(); // bottom inner (mot tunnel)
                var bo = new List<Point3d>(); // bottom outer (ytterkant nisje)
                var ti = new List<Point3d>(); // top inner
                var to = new List<Point3d>(); // top outer

                double sign = isRightSide ? +1.0 : -1.0;

                for (int i = 0; i < sampleCount; i++)
                {
                    double u = (double)i / (sampleCount - 1); // 0..1
                    double s = s0 + u * totalLen;

                    double wExtra;
                    if (s <= s1)
                    {
                        double uRamp = (s - s0) / rampLen;
                        wExtra = epsWidth + uRamp * (NisjeExtra - epsWidth);
                    }
                    else if (s <= s2)
                    {
                        wExtra = NisjeExtra;
                    }
                    else
                    {
                        double uRamp = (s3 - s) / rampLen;
                        wExtra = epsWidth + uRamp * (NisjeExtra - epsWidth);
                    }

                    if (!path.LengthParameter(s, out double t))
                    {
                        info.Add($"Fikk ikke parameter for stasjon {s:0.###} – hopper over nisje.");
                        return null;
                    }

                    if (!path.PerpendicularFrameAt(t, out Plane frame))
                    {
                        info.Add($"Fikk ikke PerpendicularFrame ved t={t:0.###} – hopper over nisje.");
                        return null;
                    }

                    Vector3d X = frame.XAxis;
                    Vector3d Y = frame.YAxis;
                    X.Unitize();
                    Y.Unitize();

                    // Som busslommer: litt inn i eksisterende vegg (InnerInset)
                    double xInner = sign * (baseOffset - InnerInset);
                    double xOuter = sign * (baseOffset + wExtra);

                    // Løft bunn litt (BottomLift) for å unngå helt koplan gulv-snitt
                    Point3d pBI = frame.Origin + X * xInner + Y * BottomLift;
                    Point3d pBO = frame.Origin + X * xOuter + Y * BottomLift;
                    Point3d pTI = frame.Origin + X * xInner + Y * (BottomLift + lommeHoyde);
                    Point3d pTO = frame.Origin + X * xOuter + Y * (BottomLift + lommeHoyde);

                    bi.Add(pBI);
                    bo.Add(pBO);
                    ti.Add(pTI);
                    to.Add(pTO);
                }

                Curve cBI = new PolylineCurve(bi);
                Curve cBO = new PolylineCurve(bo);
                Curve cTI = new PolylineCurve(ti);
                Curve cTO = new PolylineCurve(to);

                var brepParts = new List<Brep>();

                void AddLoft(Curve a, Curve b)
                {
                    var lofts = Brep.CreateFromLoft(
                        new[] { a, b },
                        Point3d.Unset,
                        Point3d.Unset,
                        LoftType.Normal,
                        false);

                    if (lofts != null && lofts.Length > 0)
                        brepParts.Add(lofts[0]);
                }

                AddLoft(cBI, cTI); // inner vegg
                AddLoft(cBO, cTO); // ytter vegg
                AddLoft(cBI, cBO); // gulv
                AddLoft(cTI, cTO); // tak

                Brep MakeCap(int idx)
                {
                    Curve[] edges = new Curve[]
                    {
                        new LineCurve(bi[idx], bo[idx]),
                        new LineCurve(bo[idx], to[idx]),
                        new LineCurve(to[idx], ti[idx]),
                        new LineCurve(ti[idx], bi[idx])
                    };
                    var caps = Brep.CreatePlanarBreps(edges, tol);
                    return (caps != null && caps.Length > 0) ? caps[0] : null;
                }

                Brep cap0 = MakeCap(0);
                Brep capN = MakeCap(bi.Count - 1);
                if (cap0 != null) brepParts.Add(cap0);
                if (capN != null) brepParts.Add(capN);

                var joined = Brep.JoinBreps(brepParts, tol);
                if (joined != null && joined.Length > 0)
                {
                    var nisje = joined[0];
                    nisje = PrepareSolid(nisje, "Nisje");
                    return nisje;
                }

                info.Add("Klarte ikke å join’e flater til en lukket nisje-brep.");
                return null;
            }

            // ---- Forbered tunnelbrep som solid (som i busslommer) ----
            Brep tunnelPrepared = PrepareSolid(tunnel.DuplicateBrep(), "Tunnel");
            if (tunnelPrepared == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Klarte ikke å forberede tunnelbrep.");
                return;
            }

            // ---- Lag nisjer på høyre og venstre side ----
            foreach (double sC in stR)
            {
                int idx = lommer.Count;
                Brep nisje = LagNisje(sC, true);
                if (nisje != null)
                {
                    nisje = PrepareSolid(nisje, $"Lomme #{idx}");
                    lommer.Add(nisje);

                    if (path.LengthParameter(sC, out double tC) &&
                        path.PerpendicularFrameAt(tC, out Plane frameC))
                    {
                        double sign = +1.0;
                        Point3d center = frameC.Origin +
                                         frameC.XAxis * (sign * (baseOffset + NisjeExtra / 2.0)) +
                                         frameC.YAxis * (BottomLift + lommeHoyde * 0.5);
                        sentre.Add(center);
                    }
                }
            }

            foreach (double sC in stL)
            {
                int idx = lommer.Count;
                Brep nisje = LagNisje(sC, false);
                if (nisje != null)
                {
                    nisje = PrepareSolid(nisje, $"Lomme #{idx}");
                    lommer.Add(nisje);

                    if (path.LengthParameter(sC, out double tC) &&
                        path.PerpendicularFrameAt(tC, out Plane frameC))
                    {
                        double sign = -1.0;
                        Point3d center = frameC.Origin +
                                         frameC.XAxis * (sign * (baseOffset + NisjeExtra / 2.0)) +
                                         frameC.YAxis * (BottomLift + lommeHoyde * 0.5);
                        sentre.Add(center);
                    }
                }
            }

            // ---- BooleanUnion pr. lomme, med flere toleranser og logging av dårlige punkter ----
            Brep tunnelOut = tunnelPrepared;
            if (lommer.Count == 0)
            {
                info.Add("Ingen nisjer generert (ingen stasjoner?).");
            }
            else
            {
                bool TryUnion(string label, Brep currentTunnel, Brep lomme, double tolTry, out Brep result)
                {
                    Point3d[] nakedEdgePts;
                    Point3d[] badIntersections;
                    Point3d[] nonManifoldEdgePts;

                    var union = Brep.CreateBooleanUnion(
                        new List<Brep> { currentTunnel, lomme },
                        tolTry,
                        manifoldOnly: false,
                        out nakedEdgePts,
                        out badIntersections,
                        out nonManifoldEdgePts);

                    int nBad = badIntersections != null ? badIntersections.Length : 0;
                    int nNaked = nakedEdgePts != null ? nakedEdgePts.Length : 0;
                    int nNonM = nonManifoldEdgePts != null ? nonManifoldEdgePts.Length : 0;

                    if (badIntersections != null) allBadPts.AddRange(badIntersections);
                    if (nakedEdgePts != null) allNakedPts.AddRange(nakedEdgePts);
                    if (nonManifoldEdgePts != null) allNonManifoldPts.AddRange(nonManifoldEdgePts);

                    if (union != null && union.Length > 0)
                    {
                        result = union[0];
                        info.Add($"{label}: BooleanUnion OK (tol={tolTry}). BadIntersectionPoints: {nBad}, NakedEdgePoints: {nNaked}, NonManifoldEdgePoints: {nNonM}.");
                        return true;
                    }

                    info.Add($"{label}: BooleanUnion FEILET (tol={tolTry}). BadIntersectionPoints: {nBad}.");
                    result = currentTunnel;
                    return false;
                }

                for (int i = 0; i < lommer.Count; i++)
                {
                    Brep lomme = lommer[i];
                    string label = $"Lomme #{i}";

                    Brep result;
                    // Først med modell-toleranse
                    if (!TryUnion(label, tunnelOut, lomme, tol, out result))
                    {
                        // Deretter strammere
                        double tolTight = tol * tolTightFactor;
                        if (!TryUnion(label + " (tight)", tunnelOut, lomme, tolTight, out result))
                        {
                            // Til slutt litt slakkere
                            double tolLoose = tol * tolLooseFactor;
                            if (!TryUnion(label + " (loose)", tunnelOut, lomme, tolLoose, out result))
                            {
                                info.Add($"{label}: BooleanUnion FEILET i alle forsøk – beholder tunnel uten denne lommen.");
                                continue;
                            }
                        }
                    }

                    tunnelOut = result;

                    // Sikre at resultatet fortsatt er brukbart
                    if (!tunnelOut.IsSolid)
                    {
                        info.Add($"{label}: Resultat etter union er ikke solid – prøver Repair.");
                        tunnelOut.Repair(tol);
                        if (tunnelOut.IsSolid)
                            info.Add($"{label}: Repair OK – union-resultat er nå solid.");
                        else
                            info.Add($"{label}: Repair kjørt men union-resultat er fortsatt ikke solid.");
                    }
                }
            }

            // ---- Output ----
            da.SetData(0, tunnelOut);
            da.SetDataList(1, lommer);
            da.SetDataList(2, sentre);
            da.SetDataList(3, info);
            da.SetDataList(4, allBadPts);
            da.SetDataList(5, allNakedPts);
            da.SetDataList(6, allNonManifoldPts);
        }

        private static readonly Dictionary<string, double> ProfileWidths = new Dictionary<string, double>
        {
            { "T9.5",  9.5 },
            { "T10.5", 10.5 },
            { "T12.5", 12.5 },
            { "T13",   13.0 },
            { "T13.5", 13.5 },
            { "T14",   14.0 },
        };


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.Stopplomme1.png"))
                {
                    return new System.Drawing.Bitmap(stream);
                }
            }
        }


        public override Guid ComponentGuid => new Guid("E5C6C069-7CC2-4E9A-9D43-0B9F9C2EDC31");
    }
}
