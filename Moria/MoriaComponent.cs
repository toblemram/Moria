using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using System.Drawing;
using System.IO;

namespace Tunnel.GH
{
    public class GH_TunnelSection_Txx_PureArcs : GH_Component
    {
        public GH_TunnelSection_Txx_PureArcs() : base(
            "TunnelProfile",
            "TunnelProfileTxx",
            "Topparc (Rh) og veggarc (Rv) for valgt T-profil (T9,5 – T14). Sweeper langs en bane.",
            "Tunnel",
            "Sections")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("ProfileType", "T", "Profiltype, f.eks. T10,5 eller T14", GH_ParamAccess.item, "T14");
            p.AddCurveParameter("Path", "P", "Kurve profilen skal sveipes langs.", GH_ParamAccess.item);
            p.AddBooleanParameter("LeftToRight", "L2R", "Tving toppbuen til å gå fra venstre (pLv) til høyre (pRv).", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddCurveParameter("Profile", "C", "PolyCurve: venstre Rv, topp Rh, høyre Rv, bunn", GH_ParamAccess.item);
            p.AddCurveParameter("Parts", "Seg", "Delkurver: venstre Rv, topp Rh, høyre Rv, bunn", GH_ParamAccess.list);
            p.AddBrepParameter("Tunnel", "B", "Sveipet Brep (tunnelvolum langs Path)", GH_ParamAccess.item);
            p.AddTextParameter("Info", "i", "Status og parametere", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string profileType = "T14";
            Curve path = null;
            bool leftToRight = true;
            da.GetData(0, ref profileType);
            da.GetData(1, ref path);
            da.GetData(2, ref leftToRight);

            var info = new List<string>();
            double tol = Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;

            // --- Profildata ---
            profileType = profileType.Replace(",", ".").ToUpper();
            if (!Profiles.TryGetValue(profileType, out double[] vals))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Ugyldig profiltype: '{profileType}'. Gyldige: {string.Join(", ", Profiles.Keys)}");
                return;
            }

            double Yv = vals[0];
            double Rv = vals[1];
            double X = vals[2];
            double Rh = vals[3];

            // --- Beregn Yh fra D = Rh - Rv ---
            double dx = X * 0.5;
            double dR = Math.Abs(Rh - Rv);
            double under = dR * dR - dx * dx;
            if (under <= -1e-9)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "(Rh - Rv)^2 < (X/2)^2. Sjekk tallene.");
                return;
            }
            double dy = Math.Sqrt(Math.Max(0.0, under));
            double Yh = Yv - dy;

            // --- Sentre ---
            Point3d CvL = new Point3d(-dx, Yv, 0);
            Point3d CvR = new Point3d(dx, Yv, 0);
            Point3d Ch = new Point3d(0, Yh, 0);

            // --- Skjæringer vegg ↔ heng (øverste) ---
            if (!CircleCircleIntersect(CvL, Rv, Ch, Rh, out Point3d pL_top, out Point3d pL_low))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Vegg/heng (venstre) skjærer ikke."); return; }
            if (!CircleCircleIntersect(CvR, Rv, Ch, Rh, out Point3d pR_top, out Point3d pR_low))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Vegg/heng (høyre) skjærer ikke."); return; }

            Point3d pLv = (pL_top.Y >= pL_low.Y) ? pL_top : pL_low;
            Point3d pRv = (pR_top.Y >= pR_low.Y) ? pR_top : pR_low;

            // --- Skjæring vegg ↔ bunn (y=0) ---
            if (!CircleLineY0Intersect(CvL, Rv, true, out Point3d bL))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Venstre vegg↔bunn skjærer ikke."); return; }
            if (!CircleLineY0Intersect(CvR, Rv, false, out Point3d bR))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Høyre vegg↔bunn skjærer ikke."); return; }

            // ---------------- TOPPARC (Rh) ----------------
            Circle cH = new Circle(new Plane(Ch, Vector3d.ZAxis), Rh);
            cH.ClosestParameter(pLv, out double tL_h);
            cH.ClosestParameter(pRv, out double tR_h);

            double dF = PositiveAngleDiff(tL_h, tR_h);
            double dB = 2.0 * Math.PI - dF;
            double midA = NormalizeAngle(tL_h + 0.5 * dF);
            double midB = NormalizeAngle(tR_h + 0.5 * dB);
            bool useForwardTop = cH.PointAt(midA).Y >= cH.PointAt(midB).Y;

            Point3d startTop = useForwardTop ? pLv : pRv;
            Point3d endTop = useForwardTop ? pRv : pLv;
            Point3d midTop = cH.PointAt(useForwardTop ? midA : midB);

            Arc arcTop = new Arc(startTop, midTop, endTop);
            Curve topArc = new ArcCurve(arcTop);
            topArc = EnsureStartsAt(topArc, leftToRight ? pLv : pRv, tol);

            // ---------------- VENSTRE Rv ----------------
            Circle cVL = new Circle(new Plane(CvL, Vector3d.ZAxis), Rv);
            cVL.ClosestParameter(pLv, out double tL1);
            cVL.ClosestParameter(bL, out double tL2);
            double midF_L = NormalizeAngle(tL1 + 0.5 * PositiveAngleDiff(tL1, tL2));
            Point3d midL = cVL.PointAt(midF_L);
            Arc arcL = new Arc(pLv, midL, bL);
            Curve leftRv = new ArcCurve(arcL);

            // ---------------- HØYRE Rv (ytterst mot +X) ----------------
            Circle cVR = new Circle(new Plane(CvR, Vector3d.ZAxis), Rv);
            cVR.ClosestParameter(pRv, out double tR1);
            cVR.ClosestParameter(bR, out double tR2);

            double dF_R = PositiveAngleDiff(tR1, tR2);
            double dB_R = 2.0 * Math.PI - dF_R;
            double midF_R = NormalizeAngle(tR1 + 0.5 * dF_R);
            double midB_R = NormalizeAngle(tR2 + 0.5 * dB_R);

            Point3d midF_R_pt = cVR.PointAt(midF_R);
            Point3d midB_R_pt = cVR.PointAt(midB_R);
            Point3d midR = (midF_R_pt.X > midB_R_pt.X) ? midF_R_pt : midB_R_pt;

            Arc arcR = new Arc(pRv, midR, bR);
            Curve rightRv = new ArcCurve(arcR);

            // ---------------- Bygg profilkurve i WorldXY ----------------
            leftRv = EnsureStartsAt(leftRv, bL, tol);
            topArc = EnsureStartsAt(topArc, pLv, tol);
            rightRv = EnsureStartsAt(rightRv, pRv, tol);

            var poly = new PolyCurve();
            poly.AppendSegment(leftRv);
            poly.AppendSegment(topArc);
            poly.AppendSegment(rightRv);

            // Bunnlinje mellom bR og bL (rett linje)
            Curve bottom = new LineCurve(new Line(bR, bL));
            poly.AppendSegment(bottom);
            poly.MakeClosed(tol);

            // ---------------- Orienter profil vinkelrett på path ----------------
            // 1) Flytt midtpunktet på bunnlinjen til origo i WorldXY
            Point3d midBottom = new Point3d(
                0.5 * (bL.X + bR.X),
                0.5 * (bL.Y + bR.Y),
                0.5 * (bL.Z + bR.Z));

            Transform toOrigin = Transform.Translation(-midBottom.X, -midBottom.Y, -midBottom.Z);
            poly.Transform(toOrigin);

            // 2) Hvis path finnes: legg profilen i plan vinkelrett på path ved start
            if (path != null)
            {
                if (path.PerpendicularFrameAt(path.Domain.T0, out Plane frame))
                {
                    // frame.Origin = path-start, X/Y/Z definert vinkelrett på path-tangent
                    Transform orient = Transform.PlaneToPlane(Plane.WorldXY, frame);
                    poly.Transform(orient);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        "Kunne ikke hente PerpendicularFrameAt fra path – profilen blir liggende i WorldXY.");
                }
            }

            // ---------------- Sweep langs path ----------------
            Brep swept = null;
            if (path != null)
            {
                var sweep = new SweepOneRail
                {
                    AngleToleranceRadians = RhinoMath.ToRadians(1.0),
                    SweepTolerance = tol
                };

                Brep[] breps = sweep.PerformSweep(path, poly);
                if (breps != null && breps.Length > 0)
                    swept = breps[0];
                else
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        "Sweep feilet – sjekk at Path er åpen/gyldig og ikke altfor kinky.");
            }

            // ---------------- Output ----------------
            da.SetData(0, poly);
            da.SetDataList(1, new List<Curve> { leftRv, topArc, rightRv, bottom });
            da.SetData(2, swept);

            info.Add($"Profil: {profileType}");
            info.Add($"Yv={Yv:0.###}, Rv={Rv:0.###}, X={X:0.###}, Rh={Rh:0.###}, Yh={Yh:0.###}");
            info.Add($"Closed={poly.IsClosed}, Sweep={(swept != null)}");
            da.SetDataList(3, info);
        }

        private static readonly Dictionary<string, double[]> Profiles = new Dictionary<string, double[]>
        {
            { "T9.5",  new[] { 1.570, 4.790, 0.450, 5.212 } },
            { "T10.5", new[] { 1.570, 4.790, 1.450, 5.950 } },
            { "T12.5", new[] { 1.570, 4.790, 3.450, 7.458 } },
            { "T13",   new[] { 1.570, 4.790, 3.950, 7.825 } },
            { "T13.5", new[] { 1.570, 4.790, 4.450, 8.053 } },
            { "T14",   new[] { 1.570, 4.790, 4.950, 8.575 } },
        };


        // =================== HJELPEFUNKSJONER ===================
        private static double NormalizeAngle(double a)
        {
            double twoPi = 2.0 * Math.PI;
            a %= twoPi;
            if (a < 0) a += twoPi;
            return a;
        }

        private static double PositiveAngleDiff(double a, double b)
        {
            a = NormalizeAngle(a);
            b = NormalizeAngle(b);
            double d = b - a;
            if (d < 0) d += 2.0 * Math.PI;
            return d;
        }

        private static bool CircleLineY0Intersect(Point3d c, double r, bool wantLeft, out Point3d p)
        {
            p = Point3d.Unset;
            double rhs = r * r - c.Y * c.Y;
            if (rhs < 0 && rhs > -1e-8) rhs = 0;
            if (rhs < 0) return false;

            double dx = Math.Sqrt(rhs);
            double x1 = c.X - dx;
            double x2 = c.X + dx;
            p = new Point3d(wantLeft ? Math.Min(x1, x2) : Math.Max(x1, x2), 0, 0);
            return true;
        }

        private static bool CircleCircleIntersect(Point3d c1, double r1, Point3d c2, double r2,
                                                  out Point3d pTop, out Point3d pLow)
        {
            pTop = pLow = Point3d.Unset;
            Vector3d d = c2 - c1;
            double D = d.Length;
            if (D < 1e-12) return false;

            double a = (r1 * r1 - r2 * r2 + D * D) / (2.0 * D);
            double hh = r1 * r1 - a * a;
            if (hh < 0 && hh > -1e-8) hh = 0;
            if (hh < 0) return false;

            Vector3d ux = d / D;
            Vector3d uy = new Vector3d(-ux.Y, ux.X, 0);

            Point3d Pm = c1 + ux * a;
            double h = Math.Sqrt(hh);

            pTop = Pm + uy * h;
            pLow = Pm - uy * h;
            return true;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.SS.png"))
                {
                    return new System.Drawing.Bitmap(stream);
                }
            }
        }


        private static Curve EnsureStartsAt(Curve c, Point3d startPt, double tol)
        {
            if (c.PointAtStart.DistanceTo(startPt) <= tol) return c;

            if (c.PointAtEnd.DistanceTo(startPt) <= tol)
            {
                Curve dup = c.DuplicateCurve();
                dup.Reverse();
                return dup;
            }

            Curve dup2 = c.DuplicateCurve();
            double dStart = c.PointAtStart.DistanceTo(startPt);
            double dEnd = c.PointAtEnd.DistanceTo(startPt);
            if (dEnd < dStart) dup2.Reverse();
            return dup2;
        }

        public override Guid ComponentGuid => new Guid("F2B2C8E7-6C2B-4A4A-9E7D-92B1AF5D7A21");


    }

}
