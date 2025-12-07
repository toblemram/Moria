using System;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Moria.TunnelGeometry
{
    /// <summary>
    /// 2D helpers for building shotcrete and insulation profiles
    /// from a tunnel T-profile cross section in WorldXY.
    /// </summary>
    public static class ShotcreteProfiles2D
    {
        /// <summary>
        /// Builds ONLY the inner shotcrete profile in 2D (open curve).
        ///
        /// - Follows the tunnel profile exactly down to Y = bottom + SkirtStartY.
        /// - At that height, on left and right, it leaves the tunnel along
        ///   a straight line with SkirtAngle and SkirtLength.
        /// - The profile is OPEN between the two skirt toes (no bottom plate).
        ///
        /// Resulting curve:
        ///   toeL -> pL -> upper (pL→pR) -> pR -> toeR
        /// </summary>
        public static bool BuildShotcreteInnerProfile2D(
            Curve tunnelProfile,
            double skirtStartY,
            double skirtAngleDeg,
            double skirtLength,
            double tol,
            out Curve shotcreteInner,
            out string error)
        {
            error = null;
            shotcreteInner = null;

            if (tunnelProfile == null || !tunnelProfile.IsClosed)
            {
                error = "Tunnel profile is null or not closed.";
                return false;
            }

            if (skirtLength <= 0.0)
            {
                error = "SkirtLength must be > 0.";
                return false;
            }

            // Bounding box to find bottom & roof
            BoundingBox bb = tunnelProfile.GetBoundingBox(true);
            double bottomY = bb.Min.Y;
            double roofY = bb.Max.Y;

            // Horizontal cut level
            double cutY = bottomY + skirtStartY;
            double epsY = tol * 10.0;

            if (cutY < bottomY + epsY) cutY = bottomY + epsY;
            if (cutY > roofY - epsY) cutY = roofY - epsY;

            // Infinite-ish cut line at cutY
            double span = bb.Diagonal.Length * 2.0;
            if (span < 10.0) span = 10.0;

            var cutLine = new LineCurve(
                new Point3d(-span, cutY, 0.0),
                new Point3d(+span, cutY, 0.0));

            var ccx = Intersection.CurveCurve(tunnelProfile, cutLine, tol, tol);
            if (ccx == null || ccx.Count < 2)
            {
                error = "Could not find two intersections between tunnel profile and SkirtStartY line.";
                return false;
            }

            // Leftmost & rightmost intersection
            Point3d pL = Point3d.Unset;
            Point3d pR = Point3d.Unset;
            double tL = 0.0, tR = 0.0;
            bool leftSet = false, rightSet = false;

            foreach (var ev in ccx)
            {
                Point3d p = ev.PointA;
                double t = ev.ParameterA;

                if (!leftSet || p.X < pL.X)
                {
                    pL = p;
                    tL = t;
                    leftSet = true;
                }

                if (!rightSet || p.X > pR.X)
                {
                    pR = p;
                    tR = t;
                    rightSet = true;
                }
            }

            if (!leftSet || !rightSet || pL.DistanceTo(pR) < tol)
            {
                error = "Could not determine distinct left/right cut points on tunnel profile.";
                return false;
            }

            // Split at tL and tR
            double[] splitParams = new double[] { tL, tR };
            Array.Sort(splitParams);

            Curve[] split = tunnelProfile.Split(splitParams);
            if (split == null || split.Length != 2)
            {
                error = "Unexpected number of pieces after splitting the tunnel profile.";
                return false;
            }

            // Pick the upper piece
            Curve upperA = split[0];
            Curve upperB = split[1];
            double yA = upperA.PointAtNormalizedLength(0.5).Y;
            double yB = upperB.PointAtNormalizedLength(0.5).Y;
            Curve upper = (yA >= yB) ? upperA : upperB;

            // Ensure upper goes from pL -> pR
            if (upper.PointAtStart.DistanceTo(pL) > upper.PointAtEnd.DistanceTo(pL))
                upper.Reverse();

            // First skirt direction (down + outward)
            double angRad = RhinoMath.ToRadians(skirtAngleDeg);
            Vector3d down = new Vector3d(0.0, -1.0, 0.0);

            Vector3d dirR = down;
            dirR.Rotate(angRad, Vector3d.ZAxis);
            dirR.Unitize();
            dirR *= skirtLength;

            Vector3d dirL = down;
            dirL.Rotate(-angRad, Vector3d.ZAxis);
            dirL.Unitize();
            dirL *= skirtLength;

            Point3d toeL = pL + dirL;
            Point3d toeR = pR + dirR;

            // Compose inner open profile
            var innerPoly = new PolyCurve();
            innerPoly.Append(new LineCurve(toeL, pL));  // left skirt
            innerPoly.Append(upper);                   // roof + upper walls
            innerPoly.Append(new LineCurve(pR, toeR)); // right skirt

            shotcreteInner = innerPoly;
            return true;
        }

        /// <summary>
        /// Builds the inner INSULATION profile in 2D (open curve),
        /// based on the inner shotcrete profile and the shotcrete thickness.
        ///
        /// Logic:
        ///  - Start from the inner shotcrete profile (toeL → pL → upper → pR → toeR).
        ///  - Offset this profile OUTWARDS by shotcreteThickness → approximate outer
        ///    surface of the concrete (same shape, bare betongtykkelse mellom).
        ///  - From the toe of this outer skirt (toeL_out, toeR_out), add an extra
        ///    second skirt defined by InsulSkirtAngle & InsulSkirtLength.
        ///
        /// Resulting curve (open):
        ///   toeL2 -> toeL1_out -> (outerShot) -> toeR1_out -> toeR2
        /// </summary>
        public static bool BuildInsulationInnerProfile2D(
            Curve shotcreteInner2D,
            double shotcreteThickness,
            double insulSkirtAngleDeg,
            double insulSkirtLength,
            double tol,
            out Curve insulationInner,
            out string error)
        {
            error = null;
            insulationInner = null;

            if (shotcreteInner2D == null || !shotcreteInner2D.IsValid)
            {
                error = "Shotcrete inner profile is null or invalid.";
                return false;
            }

            if (shotcreteThickness <= 0.0)
            {
                error = "ShotcreteThickness must be > 0.";
                return false;
            }

            if (insulSkirtLength <= 0.0)
            {
                error = "InsulSkirtLength must be > 0.";
                return false;
            }

            // 1) Offset shotcrete inner OUTWARDS by shotcreteThickness
            if (!TryOffsetOpenOutwards(
                    shotcreteInner2D,
                    shotcreteThickness,
                    tol,
                    out Curve outerShot,
                    out string offErr))
            {
                error = "Failed to offset shotcrete inner for insulation: " + offErr;
                return false;
            }

            // Start/end of outerShot are the toes of første skjørt (ytre betongoverflate)
            Point3d toeL1 = outerShot.PointAtStart;
            Point3d toeR1 = outerShot.PointAtEnd;

            // 2) Second skirt: from toeL1/toeR1 with its own angle/length
            double ang2 = RhinoMath.ToRadians(insulSkirtAngleDeg);
            Vector3d down = new Vector3d(0.0, -1.0, 0.0);

            Vector3d dirR2 = down;
            dirR2.Rotate(ang2, Vector3d.ZAxis);
            dirR2.Unitize();
            dirR2 *= insulSkirtLength;

            Vector3d dirL2 = down;
            dirL2.Rotate(-ang2, Vector3d.ZAxis);
            dirL2.Unitize();
            dirL2 *= insulSkirtLength;

            Point3d toeL2 = toeL1 + dirL2;
            Point3d toeR2 = toeR1 + dirR2;

            // 3) Compose insulation inner profile (open)
            var poly = new PolyCurve();
            poly.Append(new LineCurve(toeL2, toeL1));   // second skirt, left
            poly.Append(outerShot);                     // along outer shotcrete surface
            poly.Append(new LineCurve(toeR1, toeR2));   // second skirt, right

            insulationInner = poly;
            return true;
        }

        /// <summary>
        /// Helper: offset an open curve both +d and -d and return the one
        /// with largest bounding box diagonal (assumed to be "outwards").
        /// </summary>
        private static bool TryOffsetOpenOutwards(
            Curve inner,
            double distance,
            double tol,
            out Curve outer,
            out string error)
        {
            error = null;
            outer = null;

            if (inner == null)
            {
                error = "Inner curve is null.";
                return false;
            }

            var candidates = new System.Collections.Generic.List<Curve>();

            try
            {
                Curve[] offPlus = inner.Offset(
                    Plane.WorldXY,
                    distance,
                    tol,
                    CurveOffsetCornerStyle.Sharp);
                if (offPlus != null) candidates.AddRange(offPlus);

                Curve[] offMinus = inner.Offset(
                    Plane.WorldXY,
                    -distance,
                    tol,
                    CurveOffsetCornerStyle.Sharp);
                if (offMinus != null) candidates.AddRange(offMinus);
            }
            catch (Exception e)
            {
                error = "Offset exception: " + e.Message;
                return false;
            }

            Curve best = null;
            double bestDiag = -1.0;

            foreach (var c in candidates)
            {
                if (c == null) continue;
                BoundingBox bb = c.GetBoundingBox(true);
                double diag = bb.Diagonal.Length;
                if (diag > bestDiag)
                {
                    bestDiag = diag;
                    best = c;
                }
            }

            if (best == null)
            {
                error = "No valid offset curve found.";
                return false;
            }

            outer = best;
            return true;
        }
    }
}
