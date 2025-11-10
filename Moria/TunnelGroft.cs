using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Tunnel.GH
{
    public class GH_TunnelGroft : GH_Component
    {
        public GH_TunnelGroft() : base(
            "Tunnel.VA",
            "Grøft",
            "Lager rør og ledningsfundament + grøftform inne i tunnelen etter N500-regler.",
            "Tunnel",
            "Details")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddBrepParameter(
                "Tunnel Brep",
                "B",
                "Tunnelbrep (brukes som referanse/orientering).",
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

            p.AddIntegerParameter(
                "Antall rør",
                "N",
                "Antall rør i grøfta.",
                GH_ParamAccess.item,
                1);

            p.AddNumberParameter(
                "Rørdiametre",
                "D",
                "Liste over rørdiametre (m). DN i meter. Hvis færre enn N, brukes siste diameter om igjen.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "Dybde under bunn",
                "Depth",
                "Avstand (m) fra tunnelbunn ned til rør-senter. Positiv verdi, rør går nedover.",
                GH_ParamAccess.item,
                1.0);

            p.AddNumberParameter(
                "Relativ sideposisjon",
                "RelX",
                "-1 = venstre vegg, +1 = høyre vegg, 0 = midten. " +
                "Relativt til total tunnelbredde, så kantposisjon følger med når profiltype endres.",
                GH_ParamAccess.item,
                0.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter(
                "Rør",
                "Pipes",
                "Brep-er for alle rør.",
                GH_ParamAccess.list);

            p.AddBrepParameter(
                "Fundament",
                "Fund",
                "Ledningsfundament som Brep (lukkes om mulig).",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Grøftform",
                "Groft",
                "Grøftevolum (trapes) fra fundament opp til tunnelbunn (lukkes om mulig).",
                GH_ParamAccess.item);

            p.AddTextParameter(
                "Info",
                "i",
                "Status/diagnostikk.",
                GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Brep tunnel = null;
            Curve path = null;
            string profileType = "T14";
            int nPipes = 1;
            var diamList = new List<double>();
            double depthBelowBunn = 1.0;
            double relX = 0.0;

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
            da.GetData(3, ref nPipes);
            da.GetDataList(4, diamList);
            da.GetData(5, ref depthBelowBunn);
            da.GetData(6, ref relX);

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            var info = new List<string>();
            var pipeBreps = new List<Brep>();

            if (nPipes < 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Antall rør må være ≥ 1.");
                return;
            }

            profileType = profileType.Replace(",", ".").ToUpper();
            ProfileWidths.TryGetValue(profileType, out double Bt);
            double halfBt = Bt * 0.5;

            double curveLength = path.GetLength();
            if (curveLength <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path-kurven er for kort/ugyldig.");
                return;
            }

            if (diamList.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler rørdiametre.");
                return;
            }

            // --- klargjør diametre ---
            var diameters = new double[nPipes];
            for (int i = 0; i < nPipes; i++)
            {
                double d = diamList[Math.Min(i, diamList.Count - 1)];
                if (d <= 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Rørdiameter må være > 0.");
                    return;
                }
                diameters[i] = d;
            }

            double dMax = 0.0;
            foreach (double d in diameters)
                dMax = Math.Max(dMax, d);

            double rMax = 0.5 * dMax;

            // --- N500-regler for fundament ---
            const double minUnderkant = 0.15; // 150 mm under rørbunn
            double fundamentTykk = minUnderkant;

            // fundamentbredde: max(1.5*DN, DN + 0.20)
            double breddeFund = Math.Max(1.5 * dMax, dMax + 0.20);

            // hvis flere rør: sørg for at de får plass
            double spacing = dMax + 0.20; // senteravstand mellom rør
            if (nPipes > 1)
            {
                double totPipeWidth = (nPipes - 1) * spacing + dMax;
                breddeFund = Math.Max(breddeFund, totPipeWidth + 0.20);
            }

            // rør-senter ligger under bunn (Y negativ)
            double pipeCenterY = -Math.Abs(depthBelowBunn);

            // sjekk at det er plass til fundament under røret
            if (pipeCenterY + rMax + fundamentTykk > 0)
            {
                double minDepth = rMax + fundamentTykk + 0.05;
                pipeCenterY = -minDepth;
                info.Add($"Øker dybde slik at rør får fundament: ny dybde = {minDepth:0.###} m.");
            }

            // fundamentbunn (y mest negativ)
            double yFundBottom = pipeCenterY - rMax - fundamentTykk;
            double yFundTop = yFundBottom + fundamentTykk;   // rett under rørbunn
            double yGroftTop = 0.0;                           // tunnelbunn

            // --- beregn senter-X basert på relativ posisjon OG fundamentbredde ---
            double halfFund = 0.5 * breddeFund;

            // clamp RelX
            if (relX < -1.0) relX = -1.0;
            if (relX > 1.0) relX = 1.0;

            // ønsket posisjon ift. total bredde (ren prosent)
            double cxRequested = relX * halfBt;

            // reell maks senteroffset når fundamentet skal få plass
            double minCx = -halfBt + halfFund;
            double maxCx = halfBt - halfFund;

            double cx = 0.0;
            if (halfBt <= 0.0 || halfFund <= 0.0)
            {
                cx = 0.0;
            }
            else
            {
                if (cxRequested < minCx) cx = minCx;
                else if (cxRequested > maxCx) cx = maxCx;
                else cx = cxRequested;
            }

            // --- tverrsnitt for fundament (rektangel) i WorldXY ---
            Point3d f0 = new Point3d(cx - halfFund, yFundBottom, 0);
            Point3d f1 = new Point3d(cx + halfFund, yFundBottom, 0);
            Point3d f2 = new Point3d(cx + halfFund, yFundTop, 0);
            Point3d f3 = new Point3d(cx - halfFund, yFundTop, 0);

            var fundPoly = new Polyline(new[] { f0, f1, f2, f3, f0 });
            Curve fundSection = new PolylineCurve(fundPoly);

            // --- tverrsnitt for grøft (trapes) ---
            double depthToBunn = yGroftTop - yFundBottom; // positiv
            double breddeTop = breddeFund + 2.0 * depthToBunn;
            double halfTop = 0.5 * breddeTop;

            Point3d g0 = new Point3d(cx - halfFund, yFundBottom, 0);
            Point3d g1 = new Point3d(cx + halfFund, yFundBottom, 0);
            Point3d g2 = new Point3d(cx + halfTop, yGroftTop, 0);
            Point3d g3 = new Point3d(cx - halfTop, yGroftTop, 0);

            var groftPoly = new Polyline(new[] { g0, g1, g2, g3, g0 });
            Curve groftSection = new PolylineCurve(groftPoly);

            // --- orienter seksjoner til start-ramme og sweep langs path ---
            if (!path.PerpendicularFrameAt(path.Domain.T0, out Plane frame0))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Fikk ikke PerpendicularFrameAt på path-start.");
                return;
            }

            Transform toStart = Transform.PlaneToPlane(Plane.WorldXY, frame0);
            fundSection.Transform(toStart);
            groftSection.Transform(toStart);

            var sweep = new SweepOneRail
            {
                AngleToleranceRadians = RhinoMath.ToRadians(1.0),
                SweepTolerance = tol
            };

            Brep fundament = null;
            Brep groft = null;

            var fundBreps = sweep.PerformSweep(path, fundSection);
            if (fundBreps != null && fundBreps.Length > 0)
                fundament = fundBreps[0];
            else
                info.Add("Sweep for fundament feilet.");

            var groftBreps = sweep.PerformSweep(path, groftSection);
            if (groftBreps != null && groftBreps.Length > 0)
                groft = groftBreps[0];
            else
                info.Add("Sweep for grøftform feilet.");

            // --- CAP for å få lukket (closed) Brep på fundament og grøft ---
            if (fundament != null)
            {
                Brep cappedFund = fundament.CapPlanarHoles(tol);
                if (cappedFund != null)
                {
                    fundament = cappedFund;
                    info.Add($"Fundament lukket: IsSolid={fundament.IsSolid}");
                }
                else
                {
                    info.Add("Klarte ikke cappe fundament; forblir åpen (endeflater ikke planare?).");
                }
            }

            if (groft != null)
            {
                Brep cappedGroft = groft.CapPlanarHoles(tol);
                if (cappedGroft != null)
                {
                    groft = cappedGroft;
                    info.Add($"Grøft lukket: IsSolid={groft.IsSolid}");
                }
                else
                {
                    info.Add("Klarte ikke cappe grøft; forblir åpen (endeflater ikke planare?).");
                }
            }

            // --- rør langs path ---
            int sampleCount = 40;
            double[] tSamples = path.DivideByCount(sampleCount, true);

            for (int i = 0; i < nPipes; i++)
            {
                double d = diameters[i];
                double r = 0.5 * d;

                double xLocal = (i - (nPipes - 1) / 2.0) * spacing;
                double xTotal = cx + xLocal;

                var pts = new List<Point3d>();
                foreach (double t in tSamples)
                {
                    if (!path.PerpendicularFrameAt(t, out Plane fr))
                        continue;

                    Vector3d X = fr.XAxis;
                    Vector3d Y = fr.YAxis;
                    X.Unitize();
                    Y.Unitize();

                    Point3d p = fr.Origin + X * xTotal + Y * pipeCenterY;
                    pts.Add(p);
                }

                if (pts.Count < 2)
                {
                    info.Add("For få punkter til å lage rørbane.");
                    continue;
                }

                Curve centerCurve = Curve.CreateInterpolatedCurve(pts, 3);

                var pipeBrepsLocal = Brep.CreatePipe(
                    centerCurve,
                    r,
                    false,
                    PipeCapMode.Round,
                    true,
                    tol,
                    RhinoMath.ToRadians(1.0));

                if (pipeBrepsLocal != null && pipeBrepsLocal.Length > 0)
                    pipeBreps.Add(pipeBrepsLocal[0]);
                else
                    info.Add($"CreatePipe feilet for rør {i + 1}.");
            }

            // --- output ---
            da.SetDataList(0, pipeBreps);
            da.SetData(1, fundament);
            da.SetData(2, groft);

            info.Add($"Profil: {profileType}, Bt≈{Bt:0.###} m");
            info.Add($"RelX={relX:0.###}, ønsket X={cxRequested:0.###} m, brukt X={cx:0.###} m");
            info.Add($"Fundamentbredde={breddeFund:0.###} m, tykkelse={fundamentTykk:0.###} m");
            info.Add($"Rør-senter: dybde={Math.Abs(pipeCenterY):0.###} m under bunn");
            info.Add($"Grøft-toppbredde (ved bunn): {breddeTop:0.###} m");

            da.SetDataList(3, info);
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

        public override Guid ComponentGuid => new Guid("B8F5DD65-0E7C-4A34-AE34-9D4F3DEB8D11");
    }
}
