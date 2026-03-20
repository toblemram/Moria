using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    public class GH_TunnelTechnicalBuilding : GH_Component
    {
        public GH_TunnelTechnicalBuilding()
            : base(
                "Tunnel Technical Building",
                "TechBuilding",
                "Places rectangular technical buildings outside the tunnel along a path, trimmed against the tunnel Brep.",
                "Tunnel",
                "Tunnel geom")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddBrepParameter("Tunnel Brep", "B",
                "Tunnel Brep (may include emergency bays). Used to find the local outer edge AND to trim the building.",
                GH_ParamAccess.item);

            p.AddCurveParameter("Path", "P",
                "Centreline/path the tunnel is aligned to.",
                GH_ParamAccess.item);

            p.AddNumberParameter("Stations", "S",
                "Centre stations (m along path) for technical buildings.",
                GH_ParamAccess.list);

            p.AddBooleanParameter("RightSide", "R",
                "True = place on right side, False = left side.",
                GH_ParamAccess.item, true);

            p.AddNumberParameter("Length", "L",
                "Building length along tunnel (m).",
                GH_ParamAccess.item, 20.0);

            p.AddNumberParameter("Depth", "D",
                "Building depth out from tunnel (m).",
                GH_ParamAccess.item, 8.0);

            p.AddNumberParameter("Height", "H",
                "Building height (m).",
                GH_ParamAccess.item, 6.0);

            p.AddNumberParameter("Embed", "E",
                "How far the building penetrates into the tunnel before trimming (m). Typical 0.2–0.5.",
                GH_ParamAccess.item, 0.25);

            p.AddNumberParameter("Clearance", "C",
                "Extra offset outward from the tunnel outer edge (m).",
                GH_ParamAccess.item, 0.00);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Buildings", "Bld",
                "Trimmed building Breps.",
                GH_ParamAccess.list);

            p.AddPointParameter("Centres", "Ctr",
                "Centre points (path points) for each building station.",
                GH_ParamAccess.list);

            p.AddPlaneParameter("Frames", "F",
                "Local frames used for placement (debug).",
                GH_ParamAccess.list);

            p.AddTextParameter("Info", "i",
                "Status and diagnostics.",
                GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Brep tunnelIn = null;
            Curve path = null;
            var stations = new List<double>();
            bool rightSide = true;
            double length = 20.0;
            double depth = 8.0;
            double height = 6.0;
            double embed = 0.25;
            double clearance = 0.0;

            if (!da.GetData(0, ref tunnelIn)) return;
            if (!da.GetData(1, ref path)) return;
            da.GetDataList(2, stations);
            da.GetData(3, ref rightSide);
            da.GetData(4, ref length);
            da.GetData(5, ref depth);
            da.GetData(6, ref height);
            da.GetData(7, ref embed);
            da.GetData(8, ref clearance);

            var info = new List<string>();
            var outBreps = new List<Brep>();
            var centres = new List<Point3d>();
            var frames = new List<Plane>();

            if (tunnelIn == null || !tunnelIn.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Tunnel Brep is null or invalid.");
                return;
            }
            if (path == null || !path.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path is null or invalid.");
                return;
            }

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double pathLen = path.GetLength();
            if (pathLen <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path length is too small/invalid.");
                return;
            }

            if (length <= tol || depth <= tol || height <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Length/Depth/Height must be > 0.");
                return;
            }

            if (embed < 0) embed = 0;

            // Warn if embed is suspiciously large (like your log: embed=6)
            if (embed > 2.0)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"Embed={embed:0.###}m is very large. Typical is 0.2–0.5m.");

            // Work copy of tunnel (helps booleans)
            Brep tunnel = tunnelIn.DuplicateBrep();
            tunnel.Repair(tol);
            if (!tunnel.IsSolid)
            {
                var capped = tunnel.CapPlanarHoles(tol);
                if (capped != null) tunnel = capped;
                tunnel.Repair(tol);
            }

            info.Add($"Path length: {pathLen:0.###} m");
            info.Add($"Side: {(rightSide ? "Right" : "Left")}");
            info.Add($"Dims: L={length:0.###}, D={depth:0.###}, H={height:0.###}, embed={embed:0.###}, clearance={clearance:0.###}");

            int sideSign = rightSide ? +1 : -1;

            foreach (double s in stations)
            {
                if (s < 0.0 || s > pathLen)
                {
                    info.Add($"Station {s:0.###} out of range (0..{pathLen:0.###}) – skipped.");
                    continue;
                }

                if (!path.LengthParameter(s, out double t))
                {
                    info.Add($"Station {s:0.###}: failed to get length parameter – skipped.");
                    continue;
                }

                Point3d origin = path.PointAt(t);

                // Prefer a “upright” frame: X=right, Y≈up, Z≈tangent
                if (!TryGetUprightFrame(path, t, origin, out Plane frame))
                {
                    info.Add($"Station {s:0.###}: failed to get frame – skipped.");
                    continue;
                }

                frames.Add(frame);
                centres.Add(origin);

                // Intersect tunnel with station plane (plane normal = tangent)
                Curve[] crvs;
                Point3d[] pts;
                if (!Rhino.Geometry.Intersect.Intersection.BrepPlane(tunnel, frame, tol, out crvs, out pts) ||
                    crvs == null || crvs.Length == 0)
                {
                    info.Add($"Station {s:0.###}: no tunnel section found – skipped.");
                    continue;
                }

                Transform w2p = Transform.PlaneToPlane(frame, Plane.WorldXY);

                double bestX = rightSide ? double.NegativeInfinity : double.PositiveInfinity;
                double minY = double.PositiveInfinity;

                foreach (var c in crvs)
                {
                    if (c == null) continue;
                    Curve cc = c.DuplicateCurve();
                    cc.Transform(w2p);

                    BoundingBox bb = cc.GetBoundingBox(true);

                    // outermost X
                    double x = rightSide ? bb.Max.X : bb.Min.X;
                    if (rightSide) bestX = Math.Max(bestX, x);
                    else bestX = Math.Min(bestX, x);

                    // floor (lowest Y)
                    minY = Math.Min(minY, bb.Min.Y);
                }

                if (!double.IsFinite(bestX) || !double.IsFinite(minY))
                {
                    info.Add($"Station {s:0.###}: could not compute section extents – skipped.");
                    continue;
                }

                // Box in plane coords:
                // X lateral (right), Y up, Z along tunnel
                double xAtSurface = bestX + sideSign * clearance;
                double xInner = xAtSurface - sideSign * embed;   // into tunnel
                double xOuter = xAtSurface + sideSign * depth;   // out from tunnel

                double y0 = minY;               // sit on tunnel “floor” from section
                double y1 = y0 + height;

                double z0 = -0.5 * length;
                double z1 = +0.5 * length;

                double xmin = Math.Min(xInner, xOuter);
                double xmax = Math.Max(xInner, xOuter);

                var box = new Box(
                    Plane.WorldXY,
                    new Interval(xmin, xmax),
                    new Interval(y0, y1),
                    new Interval(z0, z1));

                Brep bld = box.ToBrep();
                if (bld == null || !bld.IsValid)
                {
                    info.Add($"Station {s:0.###}: box brep invalid – skipped.");
                    continue;
                }

                Transform p2w = Transform.PlaneToPlane(Plane.WorldXY, frame);
                bld.Transform(p2w);
                bld.Repair(tol);

                Brep trimmed = TrimBuildingAgainstTunnelRobust(bld, tunnel, tol, out string trimMsg);

                if (trimmed == null || !trimmed.IsValid)
                {
                    info.Add($"Station {s:0.###}: trim failed ({trimMsg}) – using untrimmed box.");
                    outBreps.Add(bld);
                }
                else
                {
                    trimmed.Repair(tol);
                    outBreps.Add(trimmed);
                    info.Add($"Station {s:0.###}: OK. Trim: {trimMsg}");
                }
            }

            da.SetDataList(0, outBreps);
            da.SetDataList(1, centres);
            da.SetDataList(2, frames);
            da.SetDataList(3, info);
        }

        private static bool TryGetUprightFrame(Curve path, double t, Point3d origin, out Plane frame)
        {
            frame = Plane.Unset;

            Vector3d tan = path.TangentAt(t);
            if (!tan.Unitize()) return false;

            Vector3d up = Vector3d.ZAxis;

            // right = tan x up  (gives +X as "right" in typical world coords)
            Vector3d right = Vector3d.CrossProduct(tan, up);
            if (right.Length < 1e-9)
            {
                // fallback: Rhino's own perpendicular frame
                return path.PerpendicularFrameAt(t, out frame);
            }
            right.Unitize();

            // Plane constructor will orthonormalize axes
            frame = new Plane(origin, right, up);
            return frame.IsValid;
        }

        private static Brep TrimBuildingAgainstTunnelRobust(Brep building, Brep tunnel, double tol, out string message)
        {
            message = "Trim";

            try
            {
                // 1) Try boolean difference first
                Brep[] diff = Brep.CreateBooleanDifference(building, tunnel, tol);
                if (diff != null && diff.Length > 0)
                {
                    message = "BooleanDifference OK";
                    return PickLargest(diff);
                }

                // 2) Fallback: split + keep outside pieces
                Brep[] pieces = building.Split(tunnel, tol);
                if (pieces == null || pieces.Length == 0)
                {
                    message = "Boolean failed, Split returned no pieces";
                    return null;
                }

                var keep = new List<Brep>();
                foreach (var p in pieces)
                {
                    if (p == null || !p.IsValid) continue;

                    Point3d probe = GetProbePoint(p);
                    bool inside = tunnel.IsPointInside(probe, tol, true);

                    if (!inside)
                        keep.Add(p);
                }

                if (keep.Count == 0)
                {
                    message = "Split OK, but all pieces are inside tunnel";
                    return null;
                }

                Brep[] joined = Brep.JoinBreps(keep, tol);
                if (joined != null && joined.Length > 0)
                {
                    message = $"Split+Join OK ({keep.Count} pieces)";
                    return PickLargest(joined);
                }

                message = $"Split OK (unjoined, kept {keep.Count} pieces)";
                return PickLargest(keep.ToArray());
            }
            catch (Exception ex)
            {
                message = ex.Message;
                return null;
            }
        }

        private static Point3d GetProbePoint(Brep b)
        {
            var vmp = VolumeMassProperties.Compute(b);
            if (vmp != null) return vmp.Centroid;
            return b.GetBoundingBox(true).Center;
        }

        private static Brep PickLargest(Brep[] breps)
        {
            Brep best = null;
            double bestScore = double.NegativeInfinity;

            foreach (var b in breps)
            {
                if (b == null || !b.IsValid) continue;

                double score;
                var vmp = VolumeMassProperties.Compute(b);
                if (vmp != null) score = vmp.Volume;
                else score = b.GetBoundingBox(true).Volume;

                if (score > bestScore)
                {
                    bestScore = score;
                    best = b;
                }
            }

            return best;
        }

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

        public override Guid ComponentGuid =>
            new Guid("A0E9E833-8D05-4B6A-9B7A-4B1A8B6FA6D1");
    }
}