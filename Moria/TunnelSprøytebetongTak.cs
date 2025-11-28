using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Tunnel.GH
{
    public class GH_TunnelSprøytebetongTak : GH_Component
    {
        public GH_TunnelSprøytebetongTak() : base(
            "Tunnel.Sprøytebetong (følger sving)",
            "SprøyteTak",
            "Lager et sprøytebetongtak som følger tunnellens kurvatur.",
            "Tunnel",
            "Details")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter("Path", "P", "Senterlinje/path for tunnelen.", GH_ParamAccess.item);
            p.AddTextParameter("Profiltype", "T", "T-profil: T9.5–T14", GH_ParamAccess.item, "T14");
            p.AddNumberParameter("Platehøyde", "H", "Høyde opp til toppen av betongplatene [m].", GH_ParamAccess.item, 3.0);
            p.AddNumberParameter("Tykkelse", "t", "Tykkelse på sprøytebetongen [m].", GH_ParamAccess.item, 0.25);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Tak", "B", "Sprøytebetong-tak (solid).", GH_ParamAccess.item);
            p.AddTextParameter("Info", "i", "Status og beregninger.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string profileType = "T14";
            double Hplate = 3.0, t = 0.25;

            if (!da.GetData(0, ref path) || path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler path.");
                return;
            }

            da.GetData(1, ref profileType);
            da.GetData(2, ref Hplate);
            da.GetData(3, ref t);

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            var info = new List<string>();

            profileType = profileType.Replace(",", ".").ToUpperInvariant();
            if (!Profiles.TryGetValue(profileType, out Profile prof))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Ugyldig profiltype {profileType}");
                return;
            }

            // --- Bygg profilkurver ---
            BuildProfileCurves(prof,
                out Point3d bL, out Point3d bR,
                out Point3d pLv, out Point3d pRv,
                out ArcCurve leftRv, out ArcCurve rightRv, out Curve topArc);

            double topY = Hplate + t;

            // Finn taknivå (på venstre og høyre side)
            bool okL = ArcPointAtYOnCorrectSide(leftRv, topY, true, out Point3d topL);
            bool okR = ArcPointAtYOnCorrectSide(rightRv, topY, false, out Point3d topR);
            if (!okL || !okR)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Kunne ikke finne takpunkter på profil.");
                return;
            }

            // Lag lukket tverrsnitt for takflaten
            var crvs = new List<Curve>
            {
                new LineCurve(topL, topR),
                new LineCurve(topR, pRv),
                topArc.DuplicateCurve(),
                new LineCurve(pLv, topL)
            };

            Curve[] joined = Curve.JoinCurves(crvs, tol);
            if (joined == null || joined.Length == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Kunne ikke lage lukket tverrsnitt.");
                return;
            }

            Curve section = joined[0];
            if (!section.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Tverrsnittet er ikke lukket.");
            }

            // --- Lag sweep langs path ---
            Brep[] swept = Brep.CreateFromSweep(path, section, true, tol);
            if (swept == null || swept.Length == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Sweep feilet – sjekk path og profil.");
                return;
            }

            Brep joinedBrep = (swept.Length == 1) ? swept[0] : Brep.JoinBreps(swept, tol)?[0];
            if (joinedBrep == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Kunne ikke join brep.");
                return;
            }

            Brep capped = joinedBrep.CapPlanarHoles(tol);
            if (capped == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "CapPlanarHoles feilet – åpne flater kan finnes.");
                capped = joinedBrep;
            }

            info.Add($"Sprøytebetongtak generert med tykkelse {t:0.###} m.");
            info.Add($"Starter ved platehøyde {Hplate:0.###} m, følger pathens lengde {path.GetLength():0.###} m.");
            info.Add("Takvolumet bøyes med kurvaturen og følger tunnellinjen.");

            da.SetData(0, capped);
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

            pLv = (pLtop.Y > pLlow.Y) ? pLtop : pLlow;
            pRv = (pRtop.Y < pRlow.Y) ? pRtop : pRlow;

            Plane plL = new Plane(CvL, Vector3d.ZAxis);
            Circle cL = new Circle(plL, p.Rv);
            cL.ClosestParameter(bL, out double tLb);
            cL.ClosestParameter(pLv, out double tLt);
            leftRv = new ArcCurve(cL, Math.Min(tLb, tLt), Math.Max(tLb, tLt));
            if (leftRv.PointAtEnd.Y < leftRv.PointAtStart.Y)
                leftRv.Reverse();

            Plane plR = new Plane(CvR, -Vector3d.ZAxis); plR.Flip();
            Circle cR = new Circle(plR, p.Rv);
            cR.ClosestParameter(bR, out double tRb);
            cR.ClosestParameter(pRv, out double tRt);
            rightRv = new ArcCurve(cR, Math.Min(tRb, tRt), Math.Max(tRb, tRt));
            if (rightRv.PointAtEnd.Y < rightRv.PointAtStart.Y)
                rightRv.Reverse();

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
            Point3d cand = leftSide
                ? new Point3d(Math.Min(x1, x2), yTarget, 0)
                : new Point3d(Math.Max(x1, x2), yTarget, 0);
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

        public override Guid ComponentGuid => new("A52DBE21-BF40-4A47-83B9-9E52BEFF85D4");
        protected override System.Drawing.Bitmap Icon => null;
    }
}
