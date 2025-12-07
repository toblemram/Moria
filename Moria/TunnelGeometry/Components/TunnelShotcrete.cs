using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;

using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Builds the 3D shotcrete shell along a tunnel path using a "net"
    /// profile (inner shotcrete, thickness = 0) that is swept along the
    /// path and then thickened by the user-specified Thickness.
    ///
    /// Geometry logic:
    ///  - Build inner shotcrete 2D profile:
    ///      * Follows T-profile down to SkirtStartY.
    ///      * Leaves T-profile at SkirtStartY on left/right as a straight
    ///        line with given SkirtAngle and SkirtLength.
    ///      * Profile is OPEN between the skirt toes (ingen bunnplate).
    ///  - Sweep inner profile along path (with bays → segmentert).
    ///  - For each swept surface, use Brep.CreateOffsetBrep(solid = true)
    ///    to get a constant-thickness shotcrete volume.
    /// </summary>
    public class GH_TunnelShotcrete3D : GH_Component
    {
        public GH_TunnelShotcrete3D()
            : base(
                "Tunnel Shotcrete 3D",
                "Shotcrete3D",
                "Sweeps a zero-thickness shotcrete profile along the path and thickens it to the given thickness. Supports emergency bays.",
                "Tunnel",
                "Tunnel Geom")
        {
        }

        // --------------------------------------------------------------------
        // INPUTS
        // --------------------------------------------------------------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter(
                "Path",
                "P",
                "Tunnel centreline/path.",
                GH_ParamAccess.item);

            p.AddTextParameter(
                "BaseProfileType",
                "T",
                "Base tunnel T-profile type, e.g. T5.5, T7.5, T8.5, T9.5, T10.5, T12.5, T13, T13.5, T14.",
                GH_ParamAccess.item,
                "T10.5");

            p.AddNumberParameter(
                "Thickness",
                "t",
                "Shotcrete thickness [m], applied as an offset from the inner shotcrete surface.",
                GH_ParamAccess.item,
                0.15);

            p.AddNumberParameter(
                "SkirtStartY",
                "Yk",
                "Height [m] above the tunnel bottom where the skirt starts (below this height the shotcrete leaves the T-profile).",
                GH_ParamAccess.item,
                1.5);

            p.AddNumberParameter(
                "SkirtAngle",
                "Ang",
                "Skirt angle [deg] measured from vertical downward direction, bending outwards. 0° = straight down, positive angles bend outwards.",
                GH_ParamAccess.item,
                20.0);

            p.AddNumberParameter(
                "SkirtLength",
                "L",
                "Length [m] of the straight skirt line from SkirtStartY along the given angle.",
                GH_ParamAccess.item,
                2.0);

            p.AddNumberParameter(
                "BayStationsRight",
                "SR",
                "Centre stations (m along path) for emergency bays on the right side.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "BayStationsLeft",
                "SL",
                "Centre stations (m along path) for emergency bays on the left side.",
                GH_ParamAccess.list);
        }

        // --------------------------------------------------------------------
        // OUTPUTS
        // --------------------------------------------------------------------

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter(
                "Shotcrete",
                "B",
                "Shotcrete tunnel Brep (solid, can be used for volume).",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Segment Solids",
                "Seg",
                "Shotcrete solids per segment (debug).",
                GH_ParamAccess.list);

            p.AddTextParameter(
                "Info",
                "i",
                "Status and diagnostics.",
                GH_ParamAccess.list);
        }

        // --------------------------------------------------------------------
        // INTERNAL TYPES
        // --------------------------------------------------------------------

        private class BayDef
        {
            public double SCenter;
            public int Side;      // +1 = right, -1 = left
            public double S0;     // start transition up
            public double S1;     // start full bay
            public double S2;     // end   full bay
            public double S3;     // end   transition down
        }

        private class SegmentDef
        {
            public double SStart;
            public double SEnd;

            public Curve InnerStart2D;
            public Curve InnerEnd2D;
        }

        // --------------------------------------------------------------------
        // MAIN
        // --------------------------------------------------------------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string baseProfileType = "T10.5";
            double thickness = 0.15;
            double skirtStartY = 1.5;
            double skirtAngleDeg = 20.0;
            double skirtLength = 2.0;
            var stationsRight = new List<double>();
            var stationsLeft = new List<double>();

            da.GetData(0, ref path);
            da.GetData(1, ref baseProfileType);
            da.GetData(2, ref thickness);
            da.GetData(3, ref skirtStartY);
            da.GetData(4, ref skirtAngleDeg);
            da.GetData(5, ref skirtLength);
            da.GetDataList(6, stationsRight);
            da.GetDataList(7, stationsLeft);

            var info = new List<string>();
            var segmentSolids = new List<Brep>();

            if (path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path is null.");
                return;
            }

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double pathLength = path.GetLength();

            if (pathLength <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path length is too small or invalid.");
                return;
            }

            if (thickness <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Thickness must be > 0.");
                return;
            }

            if (skirtLength <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SkirtLength must be > 0.");
                return;
            }

            skirtLength = Math.Abs(skirtLength);
            skirtAngleDeg = skirtAngleDeg % 360.0;
            baseProfileType = baseProfileType.Replace(",", ".").ToUpperInvariant();

            info.Add($"Path length: {pathLength:0.###} m");
            info.Add($"BaseProfileType: {baseProfileType}");
            info.Add($"Thickness: {thickness:0.###} m");
            info.Add($"SkirtStartY: {skirtStartY:0.###} m");
            info.Add($"SkirtAngle: {skirtAngleDeg:0.###}°");
            info.Add($"SkirtLength: {skirtLength:0.###} m");

            // ---------------------------------------------------------------
            // 1) Profile parameters
            // ---------------------------------------------------------------

            if (!ProfileType.Profiles.TryGetValue(baseProfileType, out var basePar))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Base profile type '{baseProfileType}' not found. Valid keys: {string.Join(", ", ProfileType.Profiles.Keys)}");
                return;
            }

            string bayProfileType = ResolveBayProfileType(baseProfileType);
            if (!ProfileType.Profiles.TryGetValue(bayProfileType, out var bayPar))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Bay profile type '{bayProfileType}' not found in ProfileType.Profiles.");
                return;
            }

            info.Add($"BayProfileType: {bayProfileType}");

            // ---------------------------------------------------------------
            // 2) Build tunnel profiles (WorldXY, bottom-mid at origin)
            // ---------------------------------------------------------------

            if (!BuildTunnelProfileWorld(
                    baseProfileType, basePar, tol,
                    out Curve baseTunnelWorld,
                    out string errBase))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBase);
                return;
            }

            if (!BuildTunnelProfileWorld(
                    bayProfileType, bayPar, tol,
                    out Curve bayTunnelWorld,
                    out string errBay))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBay);
                return;
            }

            // Align bay horizontally
            BoundingBox bbBase = baseTunnelWorld.GetBoundingBox(true);
            BoundingBox bbBay = bayTunnelWorld.GetBoundingBox(true);

            double xMinBase = bbBase.Min.X;
            double xMaxBase = bbBase.Max.X;
            double xMinBay = bbBay.Min.X;
            double xMaxBay = bbBay.Max.X;

            double dxRight = xMinBase - xMinBay; // align left sides (right bay)
            double dxLeft = xMaxBase - xMaxBay;  // align right sides (left bay)

            Curve bayTunnelRightWorld = bayTunnelWorld.DuplicateCurve();
            bayTunnelRightWorld.Transform(Transform.Translation(dxRight, 0, 0));

            Curve bayTunnelLeftWorld = bayTunnelWorld.DuplicateCurve();
            bayTunnelLeftWorld.Transform(Transform.Translation(dxLeft, 0, 0));

            info.Add($"Base X-span: [{xMinBase:0.###}, {xMaxBase:0.###}]");
            info.Add($"Bay  X-span: [{xMinBay:0.###}, {xMaxBay:0.###}]");
            info.Add($"Right bay dx = {dxRight:0.###}");
            info.Add($"Left  bay dx = {dxLeft:0.###}");

            // ---------------------------------------------------------------
            // 3) Build inner shotcrete 2D profiles (net profiles, thickness=0)
            // ---------------------------------------------------------------

            if (!BuildShotcreteInnerProfile2D(
                    baseTunnelWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve baseInnerWorld,
                    out string scErrBase))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Base shotcrete: " + scErrBase);
                return;
            }

            if (!BuildShotcreteInnerProfile2D(
                    bayTunnelRightWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve bayInnerRightWorld,
                    out string scErrBayR))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Right bay shotcrete: " + scErrBayR);
                return;
            }

            if (!BuildShotcreteInnerProfile2D(
                    bayTunnelLeftWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve bayInnerLeftWorld,
                    out string scErrBayL))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Left bay shotcrete: " + scErrBayL);
                return;
            }

            // ---------------------------------------------------------------
            // 4) Build bay definitions
            // ---------------------------------------------------------------

            const double rampLen = 30.0;
            const double flatLen = 30.0;

            var bays = new List<BayDef>();

            foreach (double s in stationsRight)
            {
                bays.Add(new BayDef
                {
                    SCenter = s,
                    Side = +1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                });
            }

            foreach (double s in stationsLeft)
            {
                bays.Add(new BayDef
                {
                    SCenter = s,
                    Side = -1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                });
            }

            bays.Sort((a, b) => a.SCenter.CompareTo(b.SCenter));

            foreach (var b in bays)
            {
                if (b.S0 < 0.0)
                {
                    info.Add($"Bay at s={b.SCenter:0.###} starts before path; clamping S0 to 0.");
                    b.S0 = 0.0;
                }
                if (b.S3 > pathLength)
                {
                    info.Add($"Bay at s={b.SCenter:0.###} ends after path; clamping S3 to path end.");
                    b.S3 = pathLength;
                }
            }

            info.Add($"Number of bays: {bays.Count}");

            // ---------------------------------------------------------------
            // 5) Build sweep segments (only inner profiles)
            // ---------------------------------------------------------------

            var segments = new List<SegmentDef>();
            double currentS = 0.0;
            double eps = 1e-6;

            foreach (var b in bays)
            {
                if (b.S0 > currentS + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = currentS,
                        SEnd = b.S0,
                        InnerStart2D = baseInnerWorld,
                        InnerEnd2D = baseInnerWorld
                    });
                }

                Curve innerBay = (b.Side > 0) ? bayInnerRightWorld : bayInnerLeftWorld;

                if (b.S1 > b.S0 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S0,
                        SEnd = b.S1,
                        InnerStart2D = baseInnerWorld,
                        InnerEnd2D = innerBay
                    });
                }

                if (b.S2 > b.S1 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S1,
                        SEnd = b.S2,
                        InnerStart2D = innerBay,
                        InnerEnd2D = innerBay
                    });
                }

                if (b.S3 > b.S2 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S2,
                        SEnd = b.S3,
                        InnerStart2D = innerBay,
                        InnerEnd2D = baseInnerWorld
                    });
                }

                currentS = Math.Max(currentS, b.S3);
            }

            if (currentS < pathLength - eps)
            {
                segments.Add(new SegmentDef
                {
                    SStart = currentS,
                    SEnd = pathLength,
                    InnerStart2D = baseInnerWorld,
                    InnerEnd2D = baseInnerWorld
                });
            }

            if (bays.Count == 0 && segments.Count == 0)
            {
                segments.Add(new SegmentDef
                {
                    SStart = 0.0,
                    SEnd = pathLength,
                    InnerStart2D = baseInnerWorld,
                    InnerEnd2D = baseInnerWorld
                });
            }

            info.Add($"Number of sweep segments: {segments.Count}");

            // ---------------------------------------------------------------
            // 6) Sweep inner profiles and offset to solids
            // ---------------------------------------------------------------

            var sweep = new SweepOneRail
            {
                SweepTolerance = tol,
                AngleToleranceRadians = RhinoMath.ToRadians(1.0)
            };

            int segIndex = 0;
            foreach (var seg in segments)
            {
                double segLen = seg.SEnd - seg.SStart;
                if (segLen <= tol)
                {
                    info.Add($"Segment {segIndex}: length too small, skipped.");
                    segIndex++;
                    continue;
                }

                if (!path.LengthParameter(seg.SStart, out double t0) ||
                    !path.LengthParameter(seg.SEnd, out double t1))
                {
                    info.Add($"Segment {segIndex}: could not get length parameters.");
                    segIndex++;
                    continue;
                }

                Curve rail = path.DuplicateCurve();
                rail = rail.Trim(t0, t1);
                if (rail == null || !rail.IsValid)
                {
                    info.Add($"Segment {segIndex}: failed to trim rail.");
                    segIndex++;
                    continue;
                }

                if (!path.PerpendicularFrameAt(t0, out Plane frame0) ||
                    !path.PerpendicularFrameAt(t1, out Plane frame1))
                {
                    info.Add($"Segment {segIndex}: could not get perpendicular frames.");
                    segIndex++;
                    continue;
                }

                // 3D inner profiles
                Curve c0 = seg.InnerStart2D.DuplicateCurve();
                Curve c1 = seg.InnerEnd2D.DuplicateCurve();

                Transform toStart = Transform.PlaneToPlane(Plane.WorldXY, frame0);
                Transform toEnd = Transform.PlaneToPlane(Plane.WorldXY, frame1);

                c0.Transform(toStart);
                c1.Transform(toEnd);

                Brep[] innerBreps = sweep.PerformSweep(rail, new Curve[] { c0, c1 });
                if (innerBreps == null || innerBreps.Length == 0)
                {
                    info.Add($"Segment {segIndex}: SweepOneRail failed.");
                    segIndex++;
                    continue;
                }

                Brep innerSurf = innerBreps[0];

                // Offset surface → solid with thickness
                try
                {
                    Brep[] blends, walls; // just placeholders; may be unused depending on Rhino version
                    Brep[] offset = Brep.CreateOffsetBrep(
                        innerSurf,
                        thickness,
                        true,   // solid
                        true,   // extend
                        tol,
                        out blends,
                        out walls);

                    if (offset != null && offset.Length > 0)
                    {
                        Brep solid = PrepareSolid(offset[0], $"Seg {segIndex}", tol, info);
                        segmentSolids.Add(solid);
                        info.Add($"Segment {segIndex}: sweep + offset solid OK.");
                    }
                    else
                    {
                        info.Add($"Segment {segIndex}: CreateOffsetBrep failed; using inner surface only.");
                        segmentSolids.Add(innerSurf);
                    }
                }
                catch (Exception ex)
                {
                    info.Add($"Segment {segIndex}: CreateOffsetBrep exception: {ex.Message}");
                    segmentSolids.Add(innerSurf);
                }

                segIndex++;
            }

            if (segmentSolids.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No segment solids were generated for shotcrete.");
                return;
            }

            // ---------------------------------------------------------------
            // 7) Join segments and output
            // ---------------------------------------------------------------

            Brep shotcreteTunnel;

            if (segmentSolids.Count == 1)
            {
                shotcreteTunnel = segmentSolids[0];
            }
            else
            {
                Brep[] joined = Brep.JoinBreps(segmentSolids, tol);
                if (joined != null && joined.Length > 0)
                {
                    shotcreteTunnel = joined[0];
                    info.Add($"JoinBreps: joined {segmentSolids.Count} segments into {joined.Length} Brep(s).");
                }
                else
                {
                    info.Add("JoinBreps failed – returning first segment only.");
                    shotcreteTunnel = segmentSolids[0];
                }
            }

            shotcreteTunnel = PrepareSolid(shotcreteTunnel, "ShotcreteTunnel", tol, info);

            da.SetData(0, shotcreteTunnel);
            da.SetDataList(1, segmentSolids);
            da.SetDataList(2, info);
        }

        // --------------------------------------------------------------------
        // BUILD TUNNEL PROFILE 2D
        // --------------------------------------------------------------------

        private static bool BuildTunnelProfileWorld(
            string type,
            ProfileType.ProfileParameters par,
            double tol,
            out Curve profile,
            out string error)
        {
            error = null;
            profile = null;

            PolyCurve poly;
            var dummySegments = new List<Curve>();
            var dummyDebug = new List<GeometryBase>();

            bool ok;
            if (ProfileType.IsLowRoof(type))
            {
                ok = LowRoofProfileBuilder.Build(
                    type,
                    par,
                    leftToRight: true,
                    tol: tol,
                    poly: out poly,
                    segments: out dummySegments,
                    debugGeom: out dummyDebug,
                    error: out error);
            }
            else
            {
                ok = StandardProfileBuilder.Build(
                    par,
                    leftToRight: true,
                    tol: tol,
                    poly: out poly,
                    segments: out dummySegments,
                    debugGeom: out dummyDebug,
                    error: out error);
            }

            if (!ok || poly == null)
            {
                if (string.IsNullOrWhiteSpace(error))
                    error = "Failed to build tunnel profile.";
                return false;
            }

            if (!poly.IsClosed)
            {
                error = "Tunnel profile is not closed.";
                return false;
            }

            Curve bottom = poly.SegmentCurve(poly.SegmentCount - 1);
            Point3d b0 = bottom.PointAtStart;
            Point3d b1 = bottom.PointAtEnd;
            Point3d mid = 0.5 * (b0 + b1);

            Transform toOrigin = Transform.Translation(-mid.X, -mid.Y, -mid.Z);
            poly.Transform(toOrigin);

            profile = poly;
            return true;
        }

        // --------------------------------------------------------------------
        // BUILD INNER SHOTCRETE PROFILE 2D (NETT-PROFIL, ÅPEN)
        // --------------------------------------------------------------------

        /// <summary>
        /// Builds ONLY the inner shotcrete profile in 2D (open curve).
        ///
        /// - Follows the tunnel profile exactly down to Y = bottom + SkirtStartY.
        /// - At that height, on left and right, it leaves the tunnel along
        ///   a straight line with SkirtAngle and SkirtLength.
        /// - The profile is OPEN between the two skirt toes (no bottom plate).
        /// </summary>
        private static bool BuildShotcreteInnerProfile2D(
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

            BoundingBox bb = tunnelProfile.GetBoundingBox(true);
            double bottomY = bb.Min.Y;
            double roofY = bb.Max.Y;

            double cutY = bottomY + skirtStartY;
            double epsY = tol * 10.0;

            if (cutY < bottomY + epsY) cutY = bottomY + epsY;
            if (cutY > roofY - epsY) cutY = roofY - epsY;

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

            double[] splitParams = new double[] { tL, tR };
            Array.Sort(splitParams);

            Curve[] split = tunnelProfile.Split(splitParams);
            if (split == null || split.Length != 2)
            {
                error = "Unexpected number of pieces after splitting the tunnel profile.";
                return false;
            }

            Curve upperA = split[0];
            Curve upperB = split[1];
            double yA = upperA.PointAtNormalizedLength(0.5).Y;
            double yB = upperB.PointAtNormalizedLength(0.5).Y;
            Curve upper = (yA >= yB) ? upperA : upperB;

            if (upper.PointAtStart.DistanceTo(pL) > upper.PointAtEnd.DistanceTo(pL))
                upper.Reverse();

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

            var innerPoly = new PolyCurve();
            innerPoly.Append(new LineCurve(toeL, pL)); // left skirt
            innerPoly.Append(upper);                  // roof + upper walls
            innerPoly.Append(new LineCurve(pR, toeR)); // right skirt
            // open curve
            shotcreteInner = innerPoly;

            return true;
        }

        // --------------------------------------------------------------------
        // BAY PROFILE TYPE HELPER
        // --------------------------------------------------------------------

        private static string ResolveBayProfileType(string baseProfile)
        {
            switch (baseProfile)
            {
                case "T5.5": return "T8.5";
                case "T7.5": return "T10.5";
                case "T8.5": return "T12.5";
                default: return "T14";
            }
        }

        // --------------------------------------------------------------------
        // SOLID CLEANUP
        // --------------------------------------------------------------------

        private static Brep PrepareSolid(Brep b, string label, double tol, List<string> info)
        {
            if (b == null) return null;

            if (!b.IsSolid)
            {
                Brep capped = b.CapPlanarHoles(tol);
                if (capped != null)
                {
                    b = capped;
                    info?.Add($"{label}: CapPlanarHoles applied.");
                }

                b.Repair(tol);
            }

            if (!b.IsSolid)
                info?.Add($"{label}: still not a closed solid after capping.");

            return b;
        }

        // --------------------------------------------------------------------
        // ICON & GUID
        // --------------------------------------------------------------------

        protected override Bitmap Icon
        {
            get
            {
                try
                {
                    var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                    using (var stream = assembly.GetManifestResourceStream("Moria.Resources.Shotcrete3D.png"))
                    {
                        if (stream != null)
                            return new Bitmap(stream);
                    }
                }
                catch
                {
                    // ignore
                }

                return new Bitmap(24, 24);
            }
        }

        public override Guid ComponentGuid =>
            new Guid("3E366D3F-9F1E-4A65-B3A0-6C6B9C9E8C42");
    }
}
