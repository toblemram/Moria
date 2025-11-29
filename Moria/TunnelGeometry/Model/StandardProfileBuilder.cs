using System;
using System.Collections.Generic;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Builds "standard" T-profiles (T9.5–T14)
    /// </summary>
    public static class StandardProfileBuilder
    {
        public static bool Build(
            ProfileType.ProfileParameters par,
            bool leftToRight,
            double tol,
            out PolyCurve poly,
            out List<Curve> segments,
            out List<GeometryBase> debugGeom,
            out string error)
        {
            error = null;
            poly = null;
            segments = new List<Curve>();
            debugGeom = new List<GeometryBase>();

            double Yv = par.Yv;
            double Rv = par.Rv;
            double X = par.X;
            double Rh = par.Rh;

            // ---------- Compute Yh from D = Rh - Rv (original formula) ----------
            double dx = X * 0.5;
            double dR = Math.Abs(Rh - Rv);
            double under = dR * dR - dx * dx;
            if (under <= -1e-9)
            {
                error = "(Rh - Rv)^2 < (X/2)^2. Check the profile parameters.";
                return false;
            }

            double dy = Math.Sqrt(Math.Max(0.0, under));
            double Yh = Yv - dy;

            // ---------- Centres ----------
            Point3d CvL = new Point3d(-dx, Yv, 0);
            Point3d CvR = new Point3d(dx, Yv, 0);
            Point3d Ch = new Point3d(0, Yh, 0);

            Circle cH = new Circle(new Plane(Ch, Vector3d.ZAxis), Rh);
            Circle cVL = new Circle(new Plane(CvL, Vector3d.ZAxis), Rv);
            Circle cVR = new Circle(new Plane(CvR, Vector3d.ZAxis), Rv);

            debugGeom.Add(new ArcCurve(cH));
            debugGeom.Add(new ArcCurve(cVL));
            debugGeom.Add(new ArcCurve(cVR));
            debugGeom.Add(new Point(Ch));
            debugGeom.Add(new Point(CvL));
            debugGeom.Add(new Point(CvR));

            // ---------- Intersections wall ↔ roof ----------
            if (!CircleCircleIntersect(CvL, Rv, Ch, Rh, out Point3d pL_top, out Point3d pL_low))
            {
                error = "Left wall / roof circles do not intersect.";
                return false;
            }

            if (!CircleCircleIntersect(CvR, Rv, Ch, Rh, out Point3d pR_top, out Point3d pR_low))
            {
                error = "Right wall / roof circles do not intersect.";
                return false;
            }

            Point3d pLv = (pL_top.Y >= pL_low.Y) ? pL_top : pL_low;
            Point3d pRv = (pR_top.Y >= pR_low.Y) ? pR_top : pR_low;

            debugGeom.Add(new Point(pLv));
            debugGeom.Add(new Point(pRv));

            // ---------- Intersections wall ↔ bottom (y = 0) ----------
            if (!CircleLineY0Intersect(CvL, Rv, true, out Point3d bL))
            {
                error = "Left wall / bottom do not intersect.";
                return false;
            }

            if (!CircleLineY0Intersect(CvR, Rv, false, out Point3d bR))
            {
                error = "Right wall / bottom do not intersect.";
                return false;
            }

            debugGeom.Add(new Point(bL));
            debugGeom.Add(new Point(bR));
            debugGeom.Add(new LineCurve(new Line(bL, bR)));

            // ================== ROOF ARC (Rh) ==================
            Curve topArc = MakeRoofArc(Ch, Rh, pLv, pRv, leftToRight, tol);

            // ================== LEFT WALL ARC (Rv) ==================
            Curve leftRv = MakeLeftWallArc(CvL, Rv, pLv, bL);

            // ================== RIGHT WALL ARC (Rv) ==================
            Curve rightRv = MakeRightWallArc(CvR, Rv, pRv, bR);

            leftRv = EnsureStartsAt(leftRv, bL, tol);
            topArc = EnsureStartsAt(topArc, pLv, tol);
            rightRv = EnsureStartsAt(rightRv, pRv, tol);

            // ================== Build profile ==================
            poly = new PolyCurve();
            poly.AppendSegment(leftRv);
            poly.AppendSegment(topArc);
            poly.AppendSegment(rightRv);

            Curve bottom = new LineCurve(new Line(bR, bL));
            poly.AppendSegment(bottom);
            poly.MakeClosed(tol);

            segments.Add(leftRv);
            segments.Add(topArc);
            segments.Add(rightRv);
            segments.Add(bottom);

            return true;
        }

        // --------------------------------------------------------------------
        //  Helper functions (same math as original component)
        // --------------------------------------------------------------------

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

        private static Curve MakeRoofArc(Point3d Ch, double Rh,
                                         Point3d pLv, Point3d pRv,
                                         bool leftToRight,
                                         double tol)
        {
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

            return EnsureStartsAt(topArc, leftToRight ? pLv : pRv, tol);
        }

        private static Curve MakeLeftWallArc(Point3d CvL, double Rv,
                                             Point3d pLv, Point3d bL)
        {
            Circle cVL = new Circle(new Plane(CvL, Vector3d.ZAxis), Rv);
            cVL.ClosestParameter(pLv, out double tL1);
            cVL.ClosestParameter(bL, out double tL2);

            double midF_L = NormalizeAngle(tL1 + 0.5 * PositiveAngleDiff(tL1, tL2));
            Point3d midL = cVL.PointAt(midF_L);

            Arc arcL = new Arc(pLv, midL, bL);
            return new ArcCurve(arcL);
        }

        private static Curve MakeRightWallArc(Point3d CvR, double Rv,
                                              Point3d pRv, Point3d bR)
        {
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
            return new ArcCurve(arcR);
        }
    }
}
