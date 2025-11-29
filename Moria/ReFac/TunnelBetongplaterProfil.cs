using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    public class GH_TunnelBetongplaterProfil : GH_Component
    {
        public GH_TunnelBetongplaterProfil() : base(
            "Tunnel.Betongplater",
            "PlaterProfil",
            "Betongplater som følger svingene i pathen, men forblir rette (ingen bøyning).",
            "Tunnel",
            "Details")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter("Path", "P", "Senterlinje/path for tunnelen.", GH_ParamAccess.item);
            p.AddTextParameter("Profiltype", "T", "T-profil: T9.5–T14", GH_ParamAccess.item, "T14");
            p.AddNumberParameter("Lengde L", "L", "Lengde pr. plate [m].", GH_ParamAccess.item, 2.0);
            p.AddNumberParameter("Høyde H", "H", "Høyde på platene [m].", GH_ParamAccess.item, 2.4);
            p.AddNumberParameter("Tykkelse t", "t", "Tykkelse utover fra veggbue [m].", GH_ParamAccess.item, 0.25);
            p.AddBooleanParameter("Venstre", "Lft", "Lag venstre plater.", GH_ParamAccess.item, true);
            p.AddBooleanParameter("Høyre", "Rgt", "Lag høyre plater.", GH_ParamAccess.item, true);
            p.AddNumberParameter("Spacing", "S", "Mellomrom mellom plater langs path [m].", GH_ParamAccess.item, 0.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Plater", "B", "Brep solids for betongplatene.", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Status og debug-info.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string profileType = "T14";
            double L = 2.0, H = 2.4, t = 0.25, spacing = 0.0;
            bool makeLeft = true, makeRight = true;

            if (!da.GetData(0, ref path) || path == null)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler path."); return; }

            da.GetData(1, ref profileType);
            da.GetData(2, ref L);
            da.GetData(3, ref H);
            da.GetData(4, ref t);
            da.GetData(5, ref makeLeft);
            da.GetData(6, ref makeRight);
            da.GetData(7, ref spacing);

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double pathLen = path.GetLength();

            if (pathLen <= tol)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path for kort."); return; }

            profileType = profileType.Replace(",", ".").ToUpperInvariant();
            if (!Profiles.TryGetValue(profileType, out Profile prof))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Ugyldig profiltype {profileType}"); return; }

            var info = new List<string>();
            var outBreps = new List<Brep>();

            // --- Bygg tverrsnitt ---
            BuildProfileCurves(prof, out Point3d bL, out Point3d bR, out Point3d pLv, out Point3d pRv,
                               out ArcCurve leftRv, out ArcCurve rightRv, out Curve topArc);

            double HmaxL = pLv.Y - bL.Y, HmaxR = pRv.Y - bR.Y;
            double HuseL = Math.Min(H, HmaxL), HuseR = Math.Min(H, HmaxR);

            info.Add($"Profil: {profileType}");
            info.Add($"Venstre bue L={leftRv.GetLength():0.###} m | Høyre bue L={rightRv.GetLength():0.###} m");

            Curve MakeSection(bool left)
            {
                ArcCurve rv = left ? leftRv : rightRv;
                Point3d bPt = left ? bL : bR;
                double Huse = left ? HuseL : HuseR;
                double targetY = bPt.Y + Huse;
                bool ok = ArcPointAtYOnCorrectSide(rv, targetY, left, out Point3d top);
                if (!ok) return null;

                rv.ClosestPoint(bPt, out double t0);
                rv.ClosestPoint(top, out double t1);
                Curve seg = rv.Trim(new Interval(Math.Min(t0, t1), Math.Max(t0, t1)));
                if (seg == null) seg = rv.DuplicateCurve();

                Vector3d tx = Vector3d.XAxis * (left ? -t : +t);
                Curve outer = seg.DuplicateCurve(); outer.Translate(tx);
                var bunn = new LineCurve(seg.PointAtStart, outer.PointAtStart);
                var topp = new LineCurve(seg.PointAtEnd, outer.PointAtEnd);

                var crvs = new List<Curve> { seg, topp, outer, bunn };
                Curve[] joined = Curve.JoinCurves(crvs, tol);
                if (joined == null || joined.Length == 0) return null;
                return joined[0];
            }

            Curve leftSec = makeLeft ? MakeSection(true) : null;
            Curve rightSec = makeRight ? MakeSection(false) : null;

            // --- Generer plater (rette, følger path) ---
            double pitch = L + spacing;
            int n = (int)Math.Floor(pathLen / pitch + 1e-9);

            for (int i = 0; i < n; i++)
            {
                double s0 = i * pitch;
                if (!path.LengthParameter(s0, out double t0)) continue;
                if (!path.PerpendicularFrameAt(t0, out Plane frame)) continue;

                // Bruk tangenten som lengderetning
                Vector3d dir = path.TangentAt(t0);
                dir.Unitize();
                dir *= L;

                void AddPlate(Curve sec, string side)
                {
                    if (sec == null) return;

                    // Kopi av seksjonen og plassering
                    var c = sec.DuplicateCurve();
                    c.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame));

                    // Lag en rett ekstrudering i tangentretningen
                    Surface srf = Surface.CreateExtrusion(c, dir);
                    if (srf == null)
                    {
                        info.Add($"{side}: Extrude feilet ved s0={s0:0.###}");
                        return;
                    }

                    // Gjør om til brep og lukk endene
                    Brep brep = srf.ToBrep().CapPlanarHoles(tol);
                    if (brep != null) outBreps.Add(brep);
                }

                if (makeLeft) AddPlate(leftSec, "Venstre");
                if (makeRight) AddPlate(rightSec, "Høyre");
            }

            info.Add($"Path-lengde={pathLen:0.###}, segment L={L:0.###} → {outBreps.Count} plater generert.");

            da.SetDataList(0, outBreps);
            da.SetDataList(1, info);
        }

        // ---------- Geometri fra T-profil ----------
        private static void BuildProfileCurves(Profile p,
            out Point3d bL, out Point3d bR,
            out Point3d pLv, out Point3d pRv,
            out ArcCurve leftRv, out ArcCurve rightRv,
            out Curve topArc)
        {
            double dx = 0.5 * p.X;
            Point3d CvL = new Point3d(-dx, p.Yv, 0);
            Point3d CvR = new Point3d(+dx, p.Yv, 0);
            Point3d Ch = new Point3d(0, p.Yh, 0);

            CircleLineY0Intersect(CvL, p.Rv, true, out bL);
            CircleLineY0Intersect(CvR, p.Rv, false, out bR);
            CircleCircleIntersect(CvL, p.Rv, Ch, p.Rh, out Point3d pLtop, out Point3d pLlow);
            CircleCircleIntersect(CvR, p.Rv, Ch, p.Rh, out Point3d pRtop, out Point3d pRlow);

            pLv = pLtop.Y > pLlow.Y ? pLtop : pLlow;
            pRv = pRtop.Y < pRlow.Y ? pRtop : pRlow;

            // VENSTRE
            Plane plL = new Plane(CvL, Vector3d.ZAxis);
            Circle cL = new Circle(plL, p.Rv);
            cL.ClosestParameter(bL, out double tLb);
            cL.ClosestParameter(pLv, out double tLt);
            leftRv = new ArcCurve(cL, Math.Min(tLb, tLt), Math.Max(tLb, tLt));
            if (leftRv.PointAtEnd.Y < leftRv.PointAtStart.Y) leftRv.Reverse();

            // HØYRE (speilet)
            Plane plR = new Plane(CvR, -Vector3d.ZAxis);
            plR.Flip();
            Circle cR = new Circle(plR, p.Rv);
            cR.ClosestParameter(bR, out double tRb);
            cR.ClosestParameter(pRv, out double tRt);
            rightRv = new ArcCurve(cR, Math.Min(tRb, tRt), Math.Max(tRb, tRt));
            if (rightRv.PointAtEnd.Y < rightRv.PointAtStart.Y) rightRv.Reverse();

            // TOPP
            Circle cH = new Circle(new Plane(Ch, Vector3d.ZAxis), p.Rh);
            cH.ClosestParameter(pLv, out double tH1);
            cH.ClosestParameter(pRv, out double tH2);
            topArc = new ArcCurve(cH, Math.Min(tH1, tH2), Math.Max(tH1, tH2));
        }

        private static bool ArcPointAtYOnCorrectSide(ArcCurve arc, double yTarget, bool leftSide, out Point3d pt)
        {
            pt = Point3d.Unset;
            Circle circ = new Circle(arc.Arc);
            double dy = yTarget - circ.Center.Y;
            double rhs = circ.Radius * circ.Radius - dy * dy;
            if (rhs < 0) return false;
            double dx = Math.Sqrt(rhs);
            double x1 = circ.Center.X - dx, x2 = circ.Center.X + dx;
            Point3d cand = leftSide ? new Point3d(Math.Min(x1, x2), yTarget, 0) : new Point3d(Math.Max(x1, x2), yTarget, 0);
            arc.ClosestPoint(cand, out double t);
            pt = arc.PointAt(t);
            return true;
        }

        private static bool CircleLineY0Intersect(Point3d c, double r, bool wantLeft, out Point3d p)
        {
            double rhs = r * r - c.Y * c.Y;
            double dx = Math.Sqrt(Math.Max(0, rhs));
            double x1 = c.X - dx, x2 = c.X + dx;
            p = new Point3d(wantLeft ? Math.Min(x1, x2) : Math.Max(x1, x2), 0, 0);
            return true;
        }

        private static bool CircleCircleIntersect(Point3d c1, double r1, Point3d c2, double r2,
            out Point3d pTop, out Point3d pLow)
        {
            pTop = pLow = Point3d.Unset;
            Vector3d d = c2 - c1;
            double D = d.Length;
            double a = (r1 * r1 - r2 * r2 + D * D) / (2 * D);
            double h = Math.Sqrt(Math.Max(0, r1 * r1 - a * a));
            Vector3d ux = d / D;
            Vector3d uy = new Vector3d(-ux.Y, ux.X, 0);
            Point3d Pm = c1 + ux * a;
            pTop = Pm + uy * h;
            pLow = Pm - uy * h;
            return true;
        }

        private struct Profile { public double Bt, Yv, Rv, X, Yh, Rh; }

        private static readonly Dictionary<string, Profile> Profiles = new()
        {
            { "T9.5",  new Profile{ Bt=9.5,  Yv=1.570, Rv=4.790, X=0.450,  Yh=1.213, Rh=5.212 } },
            { "T10.5", new Profile{ Bt=10.5, Yv=1.570, Rv=4.790, X=1.450,  Yh=0.664, Rh=5.950 } },
            { "T12.5", new Profile{ Bt=12.5, Yv=1.570, Rv=4.790, X=3.450,  Yh=-0.466, Rh=7.458 } },
            { "T13",   new Profile{ Bt=13.0, Yv=1.570, Rv=4.790, X=3.950,  Yh=-0.735, Rh=7.825 } },
            { "T13.5", new Profile{ Bt=13.5, Yv=1.570, Rv=4.790, X=4.450,  Yh=-0.817, Rh=8.053 } },
            { "T14",   new Profile{ Bt=14.0, Yv=1.570, Rv=4.790, X=4.950,  Yh=-1.294, Rh=8.575 } },
        };

        public override Guid ComponentGuid => new("1F283C23-5C40-47BA-9244-5FBA343E97AB");
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.Betongvegg.png"))
                {
                    return new System.Drawing.Bitmap(stream);
                }
            }
        }

    }
}
