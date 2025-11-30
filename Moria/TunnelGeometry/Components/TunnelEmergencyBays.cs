using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Builds a tunnel with integrated emergency bays by sweeping the path in segments:
    /// 
    ///   base → transition up → full bay → transition down → base
    /// 
    /// for each bay, and then joining the segment Breps.
    /// 
    /// Bay profile mapping:
    ///   T5.5  → T8.5  (bay)
    ///   T7.5  → T10.5 (bay)
    ///   T8.5  → T12.5 (bay)
    ///   T9.5+ → T14   (bay)
    /// 
    /// Right-side bay: left side of the cross-section is anchored,
    /// left-side bay: right side is anchored.
    /// 
    /// Input tunnel Brep is not used to build the geometry; the tunnel
    /// is recomputed from Path + BaseProfileType for maximum robustness.
    /// </summary>
    public class GH_TunnelEmergencyBay : GH_Component
    {
        public GH_TunnelEmergencyBay()
            : base(
                "Tunnel Emergency Bays",
                "EmergencyBays",
                "Builds a tunnel with emergency bays by sweeping path segments (no 3D booleans).",
                "Tunnel",
                "Tunnel geom")
        { }

        // --------------------------------------------------------------------
        // INPUTS
        // --------------------------------------------------------------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddBrepParameter(
                "Tunnel Brep",
                "B",
                "Base tunnel Brep (not used in this version; tunnel is rebuilt from Path and ProfileType).",
                GH_ParamAccess.item);

            p.AddCurveParameter(
                "Path",
                "P",
                "Centreline/path the tunnel is swept along.",
                GH_ParamAccess.item);

            p.AddTextParameter(
                "BaseProfileType",
                "T",
                "Base tunnel T-profile, e.g. T5.5, T7.5, T8.5, T9.5, T10.5, T12.5, T13, T13.5, T14.",
                GH_ParamAccess.item,
                "T10.5");

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
                "Tunnel With Bays",
                "B",
                "Tunnel Brep with emergency bays, built as joined sweep segments.",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Bay Segment Breps",
                "Bays",
                "Segment Breps that belong to bay regions (debug).",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "Bay Centres",
                "Ctr",
                "Centre points for each bay (on tunnel centreline at the bay station).",
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

        /// <summary>
        /// Definition of one bay in path-length coordinates.
        /// </summary>
        private class BayDef
        {
            public double SCenter; // station of bay centre
            public int Side;       // +1 = right, -1 = left
            public double S0;      // start of transition up
            public double S1;      // start of full bay
            public double S2;      // end of full bay
            public double S3;      // end of transition down
        }

        /// <summary>
        /// One rail-sweep segment along the path.
        /// </summary>
        private class SegmentDef
        {
            public double SStart;
            public double SEnd;
            public Curve ShapeStart2D;
            public Curve ShapeEnd2D;
            public bool IsBayRegion; // true if this segment belongs to a bay (transition or full)
        }

        // --------------------------------------------------------------------
        // MAIN
        // --------------------------------------------------------------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Brep inputTunnel = null; // unused geometry
            Curve path = null;
            string baseProfileType = "T10.5";
            var stationsRight = new List<double>();
            var stationsLeft = new List<double>();

            da.GetData(0, ref inputTunnel);
            da.GetData(1, ref path);
            da.GetData(2, ref baseProfileType);
            da.GetDataList(3, stationsRight);
            da.GetDataList(4, stationsLeft);

            var info = new List<string>();
            var bayCentres = new List<Point3d>();
            var baySegmentBreps = new List<Brep>();
            var badPts = new List<Point3d>();
            var nakedPts = new List<Point3d>();
            var nonManPts = new List<Point3d>();

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

            info.Add($"Path length: {pathLength:0.###} m");

            // ---------------------------------------------------------------
            // Build bay definitions (BayDef) from station lists
            // ---------------------------------------------------------------

            const double rampLen = 30.0;
            const double flatLen = 30.0;

            var bays = new List<BayDef>();

            foreach (double s in stationsRight)
            {
                var b = new BayDef
                {
                    SCenter = s,
                    Side = +1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                };
                bays.Add(b);
            }

            foreach (double s in stationsLeft)
            {
                var b = new BayDef
                {
                    SCenter = s,
                    Side = -1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                };
                bays.Add(b);
            }

            // Sort bays by centre station
            bays.Sort((a, b) => a.SCenter.CompareTo(b.SCenter));

            // Clamp and log
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

            // Bay centre points for output
            foreach (var b in bays)
            {
                if (b.SCenter < 0.0 || b.SCenter > pathLength)
                    continue;

                if (path.LengthParameter(b.SCenter, out double tC) &&
                    path.PerpendicularFrameAt(tC, out Plane frameC))
                {
                    bayCentres.Add(frameC.Origin);
                }
            }

            // ---------------------------------------------------------------
            // Resolve base and bay profile types and parameters
            // ---------------------------------------------------------------

            baseProfileType = baseProfileType.Replace(",", ".").ToUpperInvariant();

            if (!ProfileType.Profiles.TryGetValue(baseProfileType, out var baseParams))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Base profile type '{baseProfileType}' not found. Valid: {string.Join(", ", ProfileType.Profiles.Keys)}");
                return;
            }

            string bayProfileType = ResolveBayProfileType(baseProfileType);
            if (!ProfileType.Profiles.TryGetValue(bayProfileType, out var bayParams))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Bay profile type '{bayProfileType}' not found in ProfileType.Profiles.");
                return;
            }

            info.Add($"Base profile: {baseProfileType}");
            info.Add($"Bay  profile: {bayProfileType}");

            // ---------------------------------------------------------------
            // Build 2D base and bay profiles in WorldXY (bottom midpoint at origin)
            // ---------------------------------------------------------------

            if (!BuildProfileWorld(baseProfileType, baseParams, tol, out Curve baseProfileWorld, out string errBase))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBase);
                return;
            }

            if (!BuildProfileWorld(bayProfileType, bayParams, tol, out Curve bayProfileWorld, out string errBay))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBay);
                return;
            }

            BoundingBox bbBase = baseProfileWorld.GetBoundingBox(true);
            BoundingBox bbBay = bayProfileWorld.GetBoundingBox(true);

            double xMinBase = bbBase.Min.X;
            double xMaxBase = bbBase.Max.X;
            double xMinBay = bbBay.Min.X;
            double xMaxBay = bbBay.Max.X;

            // Right-side bay → left side fixed → align min X
            double dxRight = xMinBase - xMinBay;
            // Left-side bay → right side fixed → align max X
            double dxLeft = xMaxBase - xMaxBay;

            Curve bayProfileRightWorld = bayProfileWorld.DuplicateCurve();
            bayProfileRightWorld.Transform(Transform.Translation(dxRight, 0, 0));

            Curve bayProfileLeftWorld = bayProfileWorld.DuplicateCurve();
            bayProfileLeftWorld.Transform(Transform.Translation(dxLeft, 0, 0));

            info.Add($"Base X-span: [{xMinBase:0.###}, {xMaxBase:0.###}]");
            info.Add($"Bay  X-span: [{xMinBay:0.###}, {xMaxBay:0.###}]");
            info.Add($"Right bay shift dx = {dxRight:0.###} (align left side)");
            info.Add($"Left  bay shift dx = {dxLeft:0.###} (align right side)");

            // ---------------------------------------------------------------
            // Build path segments: base → transition up → bay → transition down → base
            // ---------------------------------------------------------------

            var segments = new List<SegmentDef>();
            double currentS = 0.0;
            double eps = 1e-6;

            foreach (var b in bays)
            {
                // base segment before bay
                if (b.S0 > currentS + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = currentS,
                        SEnd = b.S0,
                        ShapeStart2D = baseProfileWorld,
                        ShapeEnd2D = baseProfileWorld,
                        IsBayRegion = false
                    });
                }

                Curve bayShape2D = (b.Side > 0) ? bayProfileRightWorld : bayProfileLeftWorld;

                // transition up: base → bay
                if (b.S1 > b.S0 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S0,
                        SEnd = b.S1,
                        ShapeStart2D = baseProfileWorld,
                        ShapeEnd2D = bayShape2D,
                        IsBayRegion = true
                    });
                }

                // full bay: bay → bay
                if (b.S2 > b.S1 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S1,
                        SEnd = b.S2,
                        ShapeStart2D = bayShape2D,
                        ShapeEnd2D = bayShape2D,
                        IsBayRegion = true
                    });
                }

                // transition down: bay → base
                if (b.S3 > b.S2 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S2,
                        SEnd = b.S3,
                        ShapeStart2D = bayShape2D,
                        ShapeEnd2D = baseProfileWorld,
                        IsBayRegion = true
                    });
                }

                currentS = Math.Max(currentS, b.S3);
            }

            // final base segment after last bay
            if (currentS < pathLength - eps)
            {
                segments.Add(new SegmentDef
                {
                    SStart = currentS,
                    SEnd = pathLength,
                    ShapeStart2D = baseProfileWorld,
                    ShapeEnd2D = baseProfileWorld,
                    IsBayRegion = false
                });
            }

            // If there are no bays at all: one full base segment
            if (bays.Count == 0 && segments.Count == 0)
            {
                segments.Add(new SegmentDef
                {
                    SStart = 0.0,
                    SEnd = pathLength,
                    ShapeStart2D = baseProfileWorld,
                    ShapeEnd2D = baseProfileWorld,
                    IsBayRegion = false
                });
            }

            info.Add($"Number of sweep segments: {segments.Count}");

            // ---------------------------------------------------------------
            // Sweep each segment and join
            // ---------------------------------------------------------------

            var segmentBreps = new List<Brep>();

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

                Curve cStart2D = seg.ShapeStart2D.DuplicateCurve();
                Curve cEnd2D = seg.ShapeEnd2D.DuplicateCurve();

                Curve cStart3D = cStart2D.DuplicateCurve();
                Curve cEnd3D = cEnd2D.DuplicateCurve();

                Transform toStart = Transform.PlaneToPlane(Plane.WorldXY, frame0);
                Transform toEnd = Transform.PlaneToPlane(Plane.WorldXY, frame1);

                cStart3D.Transform(toStart);
                cEnd3D.Transform(toEnd);

                Brep[] piece = sweep.PerformSweep(rail, new Curve[] { cStart3D, cEnd3D });
                if (piece == null || piece.Length == 0)
                {
                    info.Add($"Segment {segIndex}: SweepOneRail failed.");
                    segIndex++;
                    continue;
                }

                Brep segBrep = piece[0];
                segmentBreps.Add(segBrep);
                if (seg.IsBayRegion)
                    baySegmentBreps.Add(segBrep);

                info.Add($"Segment {segIndex}: sweep OK (s={seg.SStart:0.###}→{seg.SEnd:0.###}).");
                segIndex++;
            }

            if (segmentBreps.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No segment Breps were generated.");
                return;
            }

            Brep tunnelWithBays = null;

            if (segmentBreps.Count == 1)
            {
                tunnelWithBays = segmentBreps[0];
            }
            else
            {
                Brep[] joined = Brep.JoinBreps(segmentBreps, tol);
                if (joined != null && joined.Length > 0)
                {
                    tunnelWithBays = joined[0];
                    info.Add($"JoinBreps: joined {segmentBreps.Count} segments into {joined.Length} Brep(s).");
                }
                else
                {
                    info.Add("JoinBreps failed – returning first segment only.");
                    tunnelWithBays = segmentBreps[0];
                }
            }

            // As a final polishing step, try to clean up / cap if needed
            tunnelWithBays = PrepareSolid(tunnelWithBays, "TunnelWithBays", tol, info);

            // ---------------------------------------------------------------
            // OUTPUTS
            // ---------------------------------------------------------------

            da.SetData(0, tunnelWithBays);
            da.SetDataList(1, baySegmentBreps);
            da.SetDataList(2, bayCentres);
            da.SetDataList(3, info);
            da.SetDataList(4, badPts);
            da.SetDataList(5, nakedPts);
            da.SetDataList(6, nonManPts);
        }

        // --------------------------------------------------------------------
        // Resolves which bay profile to use for a given base profile
        // --------------------------------------------------------------------

        private static string ResolveBayProfileType(string baseProfile)
        {
            switch (baseProfile)
            {
                case "T5.5": return "T8.5";
                case "T7.5": return "T10.5";
                case "T8.5": return "T12.5";
                default:
                    return "T14";
            }
        }

        // --------------------------------------------------------------------
        // Build 2D profile curve (WorldXY, bottom midpoint at origin)
        // using your Standard/LowRoof builders + ProfileType.
        // --------------------------------------------------------------------

        private bool BuildProfileWorld(
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
                return false;

            // Same convention as GH_TunnelProfile:
            // move midpoint of bottom segment to origin.
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
        // Try to make a Brep as solid/clean as possible
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

            return b;
        }

        // --------------------------------------------------------------------
        // Icon & GUID
        // --------------------------------------------------------------------

        protected override Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.Stopplomme1.png"))
                {
                    return new Bitmap(stream);
                }
            }
        }

        // Beholder samme GUID som gammel komponent slik at GH-definisjoner fortsatt virker.
        public override Guid ComponentGuid =>
            new Guid("E5C6C069-7CC2-4E9A-9D43-0B9F9C2EDC31");
    }
}
