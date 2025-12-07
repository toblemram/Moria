using System;
using System.Collections.Generic;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry
{
    /// <summary>
    /// Generates 3D rebar curves (transverse rings + longitudinal bars)
    /// based on a tunnel path and a 2D shotcrete inner profile.
    /// </summary>
    public static class TunnelRebarGenerator
    {
        public static bool GenerateRebar(
            Curve path,
            Curve innerProfile2D,
            double spacingTransverse,   // spacing around cross-section
            double spacingLongitudinal, // spacing along tunnel
            double rebarDiameter,
            double rebarCover,          // distance from inner shotcrete surface into the concrete
            double tol,
            out List<Curve> rebarCurves,
            out string error)
        {
            error = null;
            rebarCurves = new List<Curve>();

            if (path == null || !path.IsValid)
            {
                error = "Path is null or invalid.";
                return false;
            }
            if (innerProfile2D == null || !innerProfile2D.IsValid)
            {
                error = "Inner profile is null or invalid.";
                return false;
            }

            double L = path.GetLength();
            if (L <= tol)
            {
                error = "Path is too short.";
                return false;
            }

            if (spacingTransverse <= 0.0 || spacingLongitudinal <= 0.0)
            {
                error = "Rebar spacings must be > 0.";
                return false;
            }

            if (rebarDiameter <= 0.0)
            {
                error = "Rebar diameter must be > 0.";
                return false;
            }

            // -----------------------------------------------------------
            // 1) Build a 2D "rebar profile" by offsetting inner shotcrete
            //     inward by rebarCover (approximate).
            // -----------------------------------------------------------
            Curve rebar2D = innerProfile2D;
            if (rebarCover > tol)
            {
                try
                {
                    var candidates = new List<Curve>();

                    Curve[] off1 = innerProfile2D.Offset(
                        Plane.WorldXY,
                        rebarCover,
                        tol,
                        CurveOffsetCornerStyle.Sharp);
                    if (off1 != null) candidates.AddRange(off1);

                    Curve[] off2 = innerProfile2D.Offset(
                        Plane.WorldXY,
                        -rebarCover,
                        tol,
                        CurveOffsetCornerStyle.Sharp);
                    if (off2 != null) candidates.AddRange(off2);

                    Curve best = null;
                    double bestDiag = double.MaxValue;

                    // Vi velger den offsetten som havner nærmest tunnelens indre (minste bounding box)
                    foreach (var c in candidates)
                    {
                        if (c == null) continue;
                        BoundingBox bb = c.GetBoundingBox(true);
                        double diag = bb.Diagonal.Length;
                        if (diag < bestDiag)
                        {
                            bestDiag = diag;
                            best = c;
                        }
                    }

                    if (best != null)
                        rebar2D = best;
                }
                catch
                {
                    // fall-back: bruk innerProfile2D direkte
                    rebar2D = innerProfile2D;
                }
            }

            // -----------------------------------------------------------
            // 2) Transverse rings: kopier tverrsnitt langs path
            // -----------------------------------------------------------

            int nRings = (int)Math.Floor(L / spacingLongitudinal) + 1;
            for (int i = 0; i <= nRings; i++)
            {
                double s = Math.Min(i * spacingLongitudinal, L);
                if (!path.LengthParameter(s, out double t))
                    continue;

                if (!path.PerpendicularFrameAt(t, out Plane frame))
                    continue;

                Curve ring = rebar2D.DuplicateCurve();
                Transform to3D = Transform.PlaneToPlane(Plane.WorldXY, frame);
                ring.Transform(to3D);
                rebarCurves.Add(ring);
            }

            // -----------------------------------------------------------
            // 3) Longitudinal bars: punkter langs profilen → kurver langs path
            // -----------------------------------------------------------

            // Del opp rebar-profilen i punkt med gitt avstand rundt
            var pts2D = new List<Point3d>();
            double totalLen = rebar2D.GetLength();
            int nDiv = (int)Math.Floor(totalLen / spacingTransverse);
            if (nDiv < 1) nDiv = 1;

            double step = totalLen / nDiv;
            double accLen = 0.0;

            for (int i = 0; i <= nDiv; i++)
            {
                if (!rebar2D.LengthParameter(accLen, out double t2))
                    break;

                pts2D.Add(rebar2D.PointAt(t2));
                accLen += step;
            }

            // For hver 2D-posisjon lager vi en 3D-kurve som følger path
            int samplesAlong = 20; // ganske grov, men nok for visualisering
            foreach (var p2 in pts2D)
            {
                var curvePts = new List<Point3d>();

                for (int i = 0; i <= samplesAlong; i++)
                {
                    double s = L * (double)i / samplesAlong;
                    if (!path.LengthParameter(s, out double t))
                        continue;
                    if (!path.PerpendicularFrameAt(t, out Plane frame))
                        continue;

                    // Map 2D point (x,y) i WorldXY inn i dette tverrsnittet
                    double x = p2.X;
                    double y = p2.Y;
                    Point3d p3 = frame.Origin
                                 + frame.XAxis * x
                                 + frame.YAxis * y;
                    curvePts.Add(p3);
                }

                if (curvePts.Count >= 2)
                {
                    Curve bar = Curve.CreateInterpolatedCurve(curvePts, 3);
                    rebarCurves.Add(bar);
                }
            }

            return true;
        }
    }
}
