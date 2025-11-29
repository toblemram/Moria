using System;
using System.Collections.Generic;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Builds low-roof T-profiles (T5.5, T7.5, T8.5) using the special
    /// cross-construction described in the handbook N500.
    /// </summary>
    public static class LowRoofProfileBuilder
    {
        public static bool Build(
            string type,
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

            if (!ProfileType.TryGetLowRoofYh(type, out double Yh))
            {
                error = $"No Yh value defined for low-roof profile '{type}'.";
                return false;
            }

            double dx = X * 0.5;

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

            // ---------- CROSS INTERSECTION: roof ----------
            // Left roof point: intersection of Rh with RIGHT wall circle
            if (!CircleCircleIntersect(CvR, Rv, Ch, Rh, out Point3d rTop1, out Point3d rTop2))
            {
                error = "Right wall / roof circles do not intersect (low-roof).";
                return false;
            }
            Point3d pLv = PickSideThenHighest(rTop1, rTop2, wantNegative: true);

            // Right roof point: intersection of Rh with LEFT wall circle
            if (!CircleCircleIntersect(CvL, Rv, Ch, Rh, out Point3d lTop1, out Point3d lTop2))
            {
                error = "Left wall / roof circles do not intersect (low-roof).";
                return false;
            }
            Point3d pRv = PickSideThenHighest(lTop1, lTop2, wantNegative: false);

            debugGeom.Add(new Point(pLv));
            debugGeom.Add(new Point(pRv));
            debugGeom.Add(new LineCurve(new Line(Ch, pLv)));
            debugGeom.Add(new LineCurve(new Line(Ch, pRv)));

            // ---------- CROSS INTERSECTION: bottom (y=0) ----------
            // bL from RIGHT circle, inner intersection (closest to centreline)
            if (!CircleLineY0Both(CvR, Rv, out Point3d rB1, out Point3d rB2))
            {
                error = "Right wall / bottom intersections missing (low-roof).";
                return false;
            }
            Point3d bL = PickInner(rB1, rB2);   // should be x ≈ -Bt/2

            // bR from LEFT circle, inner intersection
            if (!CircleLineY0Both(CvL, Rv, out Point3d lB1, out Point3d lB2))
            {
                error = "Left wall / bottom intersections missing (low-roof).";
                return false;
            }
            Point3d bR = PickInner(lB1, lB2);   // should be x ≈ +Bt/2

            debugGeom.Add(new Point(bL));
            debugGeom.Add(new Point(bR));
            debugGeom.Add(new LineCurve(new Line(bL, bR)));

            // ---------- Build arcs ----------
            // Left wall uses RIGHT circle (cross)
            Curve leftRv = MakeInnerWallArc(cVR, bL, pLv);

            // Right wall uses LEFT circle (cross)
            Curve rightRv = MakeInnerWallArc(cVL, pRv, bR);

            // Roof arc as usual (between pLv and pRv on Rh)
            Curve topArc = MakeRoofArc(Ch, Rh, pLv, pRv, leftToRight, tol);

            leftRv = EnsureStartsAt(leftRv, bL, tol);
            topArc = EnsureStartsAt(topArc, pLv, tol);
            rightRv = EnsureStartsAt(rightRv, pRv, tol);

            // ---------- Assemble profile ----------
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
        // Helper math
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

        /// <summary>
        /// Returns both intersections of the circle with the line y=0.
        /// </summary>
        private static bool CircleLineY0Both(Point3d c, double r,
                                             out Point3d p1, out Point3d p2)
        {
            p1 = p2 = Point3d.Unset;
            double rhs = r * r - c.Y * c.Y;
            if (rhs < 0 && rhs > -1e-8) rhs = 0;
            if (rhs < 0) return false;

            double dx = Math.Sqrt(rhs);
            p1 = new Point3d(c.X - dx, 0, 0);
            p2 = new Point3d(c.X + dx, 0, 0);
            return true;
        }

        /// <summary>
        /// Picks the point with x closest to 0 (inner intersection).
        /// </summary>
        private static Point3d PickInner(Point3d a, Point3d b) =>
            (Math.Abs(a.X) <= Math.Abs(b.X)) ? a : b;

        /// <summary>
        /// Pick left or right point by x-sign, then highest Y among those.
        /// </summary>
        private static Point3d PickSideThenHighest(Point3d a, Point3d b, bool wantNegative)
        {
            bool aSide = wantNegative ? (a.X < 0) : (a.X > 0);
            bool bSide = wantNegative ? (b.X < 0) : (b.X > 0);

            if (aSide && !bSide) return a;
            if (bSide && !aSide) return b;

            // Both on same side (or none) – pick highest Y
            return (a.Y >= b.Y) ? a : b;
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

        /// <summary>
        /// Wall arc on "inner" side of the circle between start and end.
        /// Chooses the mid-arc whose midpoint is closest to the Y-axis.
        /// </summary>
        private static Curve MakeInnerWallArc(Circle circle, Point3d start, Point3d end)
        {
            circle.ClosestParameter(start, out double t1);
            circle.ClosestParameter(end, out double t2);

            double dF = PositiveAngleDiff(t1, t2);
            double dB = 2.0 * Math.PI - dF;

            double midF = NormalizeAngle(t1 + 0.5 * dF);
            double midB = NormalizeAngle(t2 + 0.5 * dB);

            Point3d midFpt = circle.PointAt(midF);
            Point3d midBpt = circle.PointAt(midB);

            Point3d mid = (Math.Abs(midFpt.X) <= Math.Abs(midBpt.X)) ? midFpt : midBpt;

            Arc arc = new Arc(start, mid, end);
            return new ArcCurve(arc);
        }
    }
}
