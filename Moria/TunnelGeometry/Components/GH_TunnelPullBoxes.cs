using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;

using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    public class GH_TunnelPullBoxes : GH_Component
    {
        private List<PullBoxInstance> _instances = new List<PullBoxInstance>();
        private PullBoxEditorWindow _editor;

        public GH_TunnelPullBoxes()
            : base(
                "Tunnel Pull Boxes",
                "PullBoxes",
                "Places electrical pull boxes along tunnel wall using Path + T-profile. Hollow box with openings on all 4 sides and thickness. Open heights are from ground (y=0) and up. Includes WPF editor.",
                "Tunnel",
                "Tunnel geom")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter("Path", "P", "Tunnel centreline/path.", GH_ParamAccess.item);

            p.AddTextParameter(
                "ProfileType", "T",
                "Tunnel T-profile, e.g. T5.5 ... T14.",
                GH_ParamAccess.item, "T10.5");

            p.AddBooleanParameter(
                "RightSide", "R",
                "If true: place at right wall. If false: left wall.",
                GH_ParamAccess.item, true);

            p.AddNumberParameter(
                "Interval", "I",
                "Spacing along centreline in metres.",
                GH_ParamAccess.item, 25.0);

            p.AddNumberParameter(
                "StartStation", "S0",
                "First station along centreline in metres.",
                GH_ParamAccess.item, 0.0);

            // Outer dims
            p.AddNumberParameter("Height", "H", "Outer height (m).", GH_ParamAccess.item, 0.50);
            p.AddNumberParameter("Width", "W", "Outer width (X) (m).", GH_ParamAccess.item, 0.58);
            p.AddNumberParameter("Length", "L", "Outer length (Z) (m).", GH_ParamAccess.item, 1.253);

            // Long-side opening (slots 8,9 stable)
            p.AddNumberParameter("LongOpenHeight", "HiL", "Opening height on LONG sides from ground (m).", GH_ParamAccess.item, 0.38);
            p.AddNumberParameter("LongOpenWidth", "WiL", "Opening width on LONG sides along tunnel (Z) (m).", GH_ParamAccess.item, 0.25);

            // Thickness (slot 10 stable)
            p.AddNumberParameter("Thickness", "Thk", "Wall/roof thickness (m).", GH_ParamAccess.item, 0.06);

            // Offsets (slots 11,12 stable)
            p.AddNumberParameter("OffsetZ", "Oz", "Vertical offset (m).", GH_ParamAccess.item, 0.0);
            p.AddNumberParameter("OffsetXY", "Oxy", "Lateral offset (m) toward/away from wall.", GH_ParamAccess.item, 0.0);

            // Optional count (slot 13 stable)
            p.AddIntegerParameter(
                "Count", "N",
                "Optional maximum number of boxes to generate. Set -1 for auto to end of path.",
                GH_ParamAccess.item, -1);

            // Short-side openings (added at end)
            p.AddNumberParameter("ShortOpenHeight", "HiS", "Opening height on SHORT sides from ground (m).", GH_ParamAccess.item, 0.38);
            p.AddNumberParameter("ShortOpenWidth", "WiS", "Opening width on SHORT sides (X) (m).", GH_ParamAccess.item, 0.25);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("PullBoxes", "B", "Pull box Breps (one per instance).", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Info/log.", GH_ParamAccess.list);

            p.AddPlaneParameter("Frames", "F", "Placement frames for each pull box.", GH_ParamAccess.list);
            p.AddNumberParameter("Stations", "S", "Station (m) for each pull box.", GH_ParamAccess.list);

            p.AddGenericParameter("PlacementInfo", "PI", "Structured info per pull box for downstream components.", GH_ParamAccess.list);
        }

        protected override void AppendAdditionalComponentMenuItems(System.Windows.Forms.ToolStripDropDown menu)
        {
            base.AppendAdditionalComponentMenuItems(menu);

            Menu_AppendSeparator(menu);

            Menu_AppendItem(menu, "Edit pull boxes…", (s, e) =>
            {
                try
                {
                    if (_editor != null)
                    {
                        _editor.Activate();
                        return;
                    }

                    _editor = new PullBoxEditorWindow(_instances);

                    _editor.ApplyRequested += (newInstances) =>
                    {
                        _instances = newInstances ?? new List<PullBoxInstance>();
                        ExpireSolution(true);
                    };

                    _editor.Closed += (ss, ee) => { _editor = null; };

                    _editor.Show(); // modeless
                }
                catch (Exception ex)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Editor error: {ex.Message}");
                }
            });
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string profileType = "T10.5";
            bool rightSide = true;

            double interval = 25.0;
            double startStation = 0.0;

            double H = 0.50, W = 0.58, L = 1.253;
            double longH = 0.38, longWz = 0.25;
            double thk = 0.06;
            double offsetZ = 0.0, offsetXY = 0.0;
            int count = -1;

            double shortH = 0.38, shortWx = 0.25;

            da.GetData(0, ref path);
            da.GetData(1, ref profileType);
            da.GetData(2, ref rightSide);
            da.GetData(3, ref interval);
            da.GetData(4, ref startStation);

            da.GetData(5, ref H);
            da.GetData(6, ref W);
            da.GetData(7, ref L);

            da.GetData(8, ref longH);
            da.GetData(9, ref longWz);
            da.GetData(10, ref thk);

            da.GetData(11, ref offsetZ);
            da.GetData(12, ref offsetXY);

            da.GetData(13, ref count);

            da.GetData(14, ref shortH);
            da.GetData(15, ref shortWx);

            var info = new List<string>();

            if (path == null || !path.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path is null or invalid.");
                return;
            }

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double pathLength = path.GetLength();
            if (pathLength <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path length too small/invalid.");
                return;
            }

            profileType = profileType.Replace(",", ".").ToUpperInvariant();

            if (!ProfileType.Profiles.TryGetValue(profileType, out var par))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Invalid profile type '{profileType}'. Valid: {string.Join(", ", ProfileType.Profiles.Keys)}");
                return;
            }

            if (!BuildProfileWorld(profileType, par, tol, out Curve profileWorld, out string err))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, err);
                return;
            }

            BoundingBox bb = profileWorld.GetBoundingBox(true);
            double xWall = rightSide ? bb.Max.X : bb.Min.X;
            info.Add($"Profile: {profileType}, xWall={xWall:0.###} (profile WorldXY space)");

            EnsureInstances(pathLength, startStation, interval, count,
                H, W, L,
                longWz, longH,
                shortWx, shortH,
                thk,
                offsetZ, offsetXY,
                info);

            UpdateStationsFromInputs(startStation, interval);
            ApplyDefaultsToInstances(
                H, W, L,
                longWz, longH,
                shortWx, shortH,
                thk,
                offsetZ, offsetXY);

            var breps = new List<Brep>();
            var frames = new List<Plane>();
            var stations = new List<double>();
            var placementInfo = new List<PullBoxPlacementInfo>();

            foreach (var inst in _instances.OrderBy(x => x.EffectiveStation))
            {
                double s = inst.EffectiveStation;
                if (s < 0 || s > pathLength)
                    continue;

                if (!path.LengthParameter(s, out double t))
                    continue;

                if (!TryMakeStationFrame(path, t, rightSide, out Plane stationFrame))
                    continue;

                Plane crossPlane = stationFrame;

                double sideSign = rightSide ? +1.0 : -1.0;
                crossPlane.Origin += crossPlane.XAxis * (inst.OffsetXY * sideSign);
                crossPlane.Origin += crossPlane.YAxis * inst.OffsetZ;

                Brep b = BuildPullBoxBrep(
                    crossPlane,
                    xWall,
                    rightSide,
                    inst.Length,
                    inst.Width,
                    inst.Height,
                    inst.LongOpenWidth,
                    inst.LongOpenHeight,
                    inst.ShortOpenWidth,
                    inst.ShortOpenHeight,
                    inst.Thickness,
                    tol,
                    out string buildMsg);

                if (!string.IsNullOrWhiteSpace(buildMsg))
                    info.Add(buildMsg);

                if (b != null && b.IsValid)
                {
                    breps.Add(b);
                    frames.Add(crossPlane);
                    stations.Add(s);

                    placementInfo.Add(new PullBoxPlacementInfo
                    {
                        Id = inst.Id,
                        Station = s,
                        Frame = crossPlane,

                        Length = inst.Length,
                        Width = inst.Width,
                        Height = inst.Height,

                        LongOpenWidth = inst.LongOpenWidth,
                        LongOpenHeight = inst.LongOpenHeight,

                        ShortOpenWidth = inst.ShortOpenWidth,
                        ShortOpenHeight = inst.ShortOpenHeight,

                        Thickness = inst.Thickness,

                        OffsetZ = inst.OffsetZ,
                        OffsetXY = inst.OffsetXY,

                        RightSide = rightSide
                    });
                }
            }

            info.Add($"Instances: {_instances.Count}, OutputBreps: {breps.Count}");

            da.SetDataList(0, breps);
            da.SetDataList(1, info);
            da.SetDataList(2, frames);
            da.SetDataList(3, stations);
            da.SetDataList(4, placementInfo.Select(pi => new GH_ObjectWrapper(pi)));
        }

        private void EnsureInstances(
            double pathLength,
            double startStation,
            double interval,
            int count,
            double H, double W, double L,
            double longWz, double longH,
            double shortWx, double shortH,
            double thk,
            double offsetZ, double offsetXY,
            List<string> info)
        {
            if (_instances != null && _instances.Count > 0)
                return;

            _instances = new List<PullBoxInstance>();

            if (interval <= 1e-9) interval = 25.0;

            int maxN;
            if (count > 0) maxN = count;
            else maxN = (int)Math.Floor(Math.Max(0.0, (pathLength - startStation)) / interval) + 1;

            for (int i = 0; i < maxN; i++)
            {
                double s = startStation + i * interval;
                if (s > pathLength) break;

                _instances.Add(new PullBoxInstance
                {
                    Station = s,
                    FollowInterval = true,
                    UseDefaults = true,

                    Height = H,
                    Width = W,
                    Length = L,

                    LongOpenWidth = longWz,
                    LongOpenHeight = longH,

                    ShortOpenWidth = shortWx,
                    ShortOpenHeight = shortH,

                    Thickness = thk,

                    OffsetZ = offsetZ,
                    OffsetXY = offsetXY
                });
            }

            info.Add($"Generated default instances: {_instances.Count} (start={startStation:0.###}, interval={interval:0.###})");
        }

        private void UpdateStationsFromInputs(double startStation, double interval)
        {
            if (_instances == null || _instances.Count == 0) return;
            if (interval <= 1e-9) return;

            var ordered = _instances.OrderBy(x => x.EffectiveStation).ToList();

            double nextS = startStation;

            foreach (var inst in ordered)
            {
                if (inst.FollowInterval)
                    inst.Station = nextS;

                nextS = inst.EffectiveStation + interval;
            }

            _instances = ordered;
        }

        private void ApplyDefaultsToInstances(
            double H, double W, double L,
            double longWz, double longH,
            double shortWx, double shortH,
            double thk,
            double offsetZ, double offsetXY)
        {
            if (_instances == null) return;

            foreach (var inst in _instances)
            {
                if (!inst.UseDefaults) continue;

                inst.Height = H;
                inst.Width = W;
                inst.Length = L;

                inst.LongOpenWidth = longWz;
                inst.LongOpenHeight = longH;

                inst.ShortOpenWidth = shortWx;
                inst.ShortOpenHeight = shortH;

                inst.Thickness = thk;

                inst.OffsetZ = offsetZ;
                inst.OffsetXY = offsetXY;
            }
        }

        private static bool TryMakeStationFrame(Curve path, double t, bool rightSide, out Plane frame)
        {
            frame = Plane.Unset;

            Point3d pt = path.PointAt(t);
            Vector3d tan = path.TangentAt(t);
            if (!tan.Unitize()) return false;

            Vector3d up = Vector3d.ZAxis;

            if (Math.Abs(Vector3d.Multiply(tan, up)) > 0.999)
            {
                if (path.PerpendicularFrameAt(t, out Plane p))
                {
                    frame = p;
                    return true;
                }
                return false;
            }

            Vector3d side = Vector3d.CrossProduct(tan, up);
            if (!side.Unitize()) return false;

            if (!rightSide) side *= -1.0;

            frame = new Plane(pt, side, up);
            return true;
        }

        /// <summary>
        /// Hollow box: outer - cavity.
        /// Openings on ALL 4 vertical faces:
        ///  - BOTH long sides (X faces): LongOpenWidth along Z, LongOpenHeight from y=0
        ///  - BOTH short sides (Z faces): ShortOpenWidth along X, ShortOpenHeight from y=0
        /// Local coords in crossPlane: X=side, Y=up, Z=along tunnel.
        /// </summary>
        /// <summary>
        /// Hollow box: outer - cavity.
        /// Openings on ALL 4 vertical faces:
        ///  - BOTH long sides (X faces)
        ///  - BOTH short sides (Z faces)
        /// Open heights are measured from ground (y=0).
        /// </summary>
        private static Brep BuildPullBoxBrep(
            Plane crossPlane,
            double xWallProfileSpace,
            bool rightSide,
            double lengthAlongTunnel,
            double widthFromWallTowardCentre,
            double height,
            double longOpenWidthZ,
            double longOpenHeightFromGround,
            double shortOpenWidthX,
            double shortOpenHeightFromGround,
            double thickness,
            double tol,
            out string info)
        {
            info = null;
            thickness = Math.Max(0.0, thickness);

            double xWall = xWallProfileSpace;

            // Outer extents
            double x0 = rightSide ? (xWall - widthFromWallTowardCentre) : xWall;
            double x1 = rightSide ? xWall : (xWall + widthFromWallTowardCentre);

            Interval ix = new Interval(Math.Min(x0, x1), Math.Max(x0, x1));
            Interval iy = new Interval(0, height);
            Interval iz = new Interval(-lengthAlongTunnel * 0.5, lengthAlongTunnel * 0.5);

            Brep outer = new Box(crossPlane, ix, iy, iz).ToBrep();
            if (outer == null || !outer.IsValid)
            {
                info = "Outer box creation failed.";
                return null;
            }

            // -------------------------
            // Hollow cavity
            // -------------------------
            Brep shell = outer;

            double xi0 = ix.T0 + thickness;
            double xi1 = ix.T1 - thickness;
            double yi0 = iy.T0 + thickness;
            double yi1 = iy.T1 - thickness;
            double zi0 = iz.T0 + thickness;
            double zi1 = iz.T1 - thickness;

            if (xi1 > xi0 + tol && yi1 > yi0 + tol && zi1 > zi0 + tol)
            {
                Brep cavity = new Box(crossPlane,
                    new Interval(xi0, xi1),
                    new Interval(yi0, yi1),
                    new Interval(zi0, zi1)).ToBrep();

                if (cavity != null && cavity.IsValid)
                {
                    Brep[] diff = Brep.CreateBooleanDifference(outer, cavity, tol);
                    if (diff != null && diff.Length > 0)
                        shell = diff[0];
                    else
                        info = "Hollow boolean failed.";
                }
            }
            else
            {
                info = "Thickness too large -> cavity invalid.";
            }

            double eps = Math.Max(tol * 10.0, 1e-4);

            // ============================================================
            // LONG SIDE OPENINGS (both X faces)
            // ============================================================
            if (longOpenWidthZ > tol && longOpenHeightFromGround > tol && thickness > tol)
            {
                double oy0 = 0.0;
                double oy1 = Math.Min(height, longOpenHeightFromGround);

                double oz0 = -longOpenWidthZ * 0.5;
                double oz1 = +longOpenWidthZ * 0.5;

                oz0 = Math.Max(iz.T0 + eps, oz0);
                oz1 = Math.Min(iz.T1 - eps, oz1);

                if (oy1 > oy0 + tol && oz1 > oz0 + tol)
                {
                    // X min face
                    {
                        double faceX = ix.T0;
                        double cx0 = faceX - eps;
                        double cx1 = faceX + thickness + eps;

                        Brep cut = new Box(crossPlane,
                            new Interval(Math.Min(cx0, cx1), Math.Max(cx0, cx1)),
                            new Interval(oy0, oy1),
                            new Interval(oz0, oz1)).ToBrep();

                        if (cut != null && cut.IsValid)
                        {
                            Brep[] d = Brep.CreateBooleanDifference(shell, cut, tol);
                            if (d != null && d.Length > 0)
                                shell = d[0];
                        }
                    }

                    // X max face
                    {
                        double faceX = ix.T1;
                        double cx0 = faceX - thickness - eps;
                        double cx1 = faceX + eps;

                        Brep cut = new Box(crossPlane,
                            new Interval(Math.Min(cx0, cx1), Math.Max(cx0, cx1)),
                            new Interval(oy0, oy1),
                            new Interval(oz0, oz1)).ToBrep();

                        if (cut != null && cut.IsValid)
                        {
                            Brep[] d = Brep.CreateBooleanDifference(shell, cut, tol);
                            if (d != null && d.Length > 0)
                                shell = d[0];
                        }
                    }
                }
            }

            // ============================================================
            // SHORT SIDE OPENINGS (both Z faces)
            // ============================================================
            if (shortOpenWidthX > tol && shortOpenHeightFromGround > tol && thickness > tol)
            {
                double cy0 = 0.0;
                double cy1 = Math.Min(height, shortOpenHeightFromGround);

                double xMid = 0.5 * (ix.T0 + ix.T1);
                double cx0 = xMid - shortOpenWidthX * 0.5;
                double cx1 = xMid + shortOpenWidthX * 0.5;

                cx0 = Math.Max(ix.T0 + eps, cx0);
                cx1 = Math.Min(ix.T1 - eps, cx1);

                if (cy1 > cy0 + tol && cx1 > cx0 + tol)
                {
                    // Z min face
                    {
                        double faceZ = iz.T0;
                        double cz0 = faceZ - eps;
                        double cz1 = faceZ + thickness + eps;

                        Brep cut = new Box(crossPlane,
                            new Interval(cx0, cx1),
                            new Interval(cy0, cy1),
                            new Interval(Math.Min(cz0, cz1), Math.Max(cz0, cz1))).ToBrep();

                        if (cut != null && cut.IsValid)
                        {
                            Brep[] d = Brep.CreateBooleanDifference(shell, cut, tol);
                            if (d != null && d.Length > 0)
                                shell = d[0];
                        }
                    }

                    // Z max face
                    {
                        double faceZ = iz.T1;
                        double cz0 = faceZ - thickness - eps;
                        double cz1 = faceZ + eps;

                        Brep cut = new Box(crossPlane,
                            new Interval(cx0, cx1),
                            new Interval(cy0, cy1),
                            new Interval(Math.Min(cz0, cz1), Math.Max(cz0, cz1))).ToBrep();

                        if (cut != null && cut.IsValid)
                        {
                            Brep[] d = Brep.CreateBooleanDifference(shell, cut, tol);
                            if (d != null && d.Length > 0)
                                shell = d[0];
                        }
                    }
                }
            }

            return shell;
        }

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
                    type, par, leftToRight: true, tol: tol,
                    poly: out poly, segments: out dummySegments, debugGeom: out dummyDebug, error: out error);
            }
            else
            {
                ok = StandardProfileBuilder.Build(
                    par, leftToRight: true, tol: tol,
                    poly: out poly, segments: out dummySegments, debugGeom: out dummyDebug, error: out error);
            }

            if (!ok || poly == null)
                return false;

            Curve bottom = poly.SegmentCurve(poly.SegmentCount - 1);
            Point3d b0 = bottom.PointAtStart;
            Point3d b1 = bottom.PointAtEnd;
            Point3d mid = 0.5 * (b0 + b1);

            Transform toOrigin = Transform.Translation(-mid.X, -mid.Y, -mid.Z);
            poly.Transform(toOrigin);

            profile = poly;
            return true;
        }

        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetInt32("InstanceCount", _instances?.Count ?? 0);

            if (_instances != null)
            {
                for (int i = 0; i < _instances.Count; i++)
                {
                    var it = _instances[i];
                    var chunk = writer.CreateChunk("Instance", i);

                    chunk.SetString("Id", it.Id.ToString());

                    chunk.SetDouble("Station", it.Station);
                    chunk.SetDouble("StationOffset", it.StationOffset);
                    chunk.SetBoolean("FollowInterval", it.FollowInterval);
                    chunk.SetBoolean("UseDefaults", it.UseDefaults);

                    chunk.SetDouble("Length", it.Length);
                    chunk.SetDouble("Width", it.Width);
                    chunk.SetDouble("Height", it.Height);

                    chunk.SetDouble("LongOpenWidth", it.LongOpenWidth);
                    chunk.SetDouble("LongOpenHeight", it.LongOpenHeight);

                    chunk.SetDouble("ShortOpenWidth", it.ShortOpenWidth);
                    chunk.SetDouble("ShortOpenHeight", it.ShortOpenHeight);

                    chunk.SetDouble("Thickness", it.Thickness);

                    chunk.SetDouble("OffsetZ", it.OffsetZ);
                    chunk.SetDouble("OffsetXY", it.OffsetXY);
                }
            }

            return base.Write(writer);
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            int n = reader.GetInt32("InstanceCount");
            _instances = new List<PullBoxInstance>();

            for (int i = 0; i < n; i++)
            {
                var chunk = reader.FindChunk("Instance", i);
                if (chunk == null) continue;

                var it = new PullBoxInstance();

                if (Guid.TryParse(chunk.GetString("Id"), out Guid gid))
                    it.Id = gid;

                it.Station = chunk.GetDouble("Station");
                it.StationOffset = chunk.GetDouble("StationOffset");
                it.FollowInterval = chunk.GetBoolean("FollowInterval");
                it.UseDefaults = chunk.GetBoolean("UseDefaults");

                it.Length = chunk.GetDouble("Length");
                it.Width = chunk.GetDouble("Width");
                it.Height = chunk.GetDouble("Height");

                double low = chunk.GetDouble("LongOpenWidth");
                double loh = chunk.GetDouble("LongOpenHeight");

                // Backwards compatibility: old InnerWidth/InnerHeight -> treat as LongOpen*
                if (low <= 0) low = chunk.GetDouble("InnerWidth");
                if (loh <= 0) loh = chunk.GetDouble("InnerHeight");

                it.LongOpenWidth = low > 0 ? low : it.LongOpenWidth;
                it.LongOpenHeight = loh > 0 ? loh : it.LongOpenHeight;

                double sow = chunk.GetDouble("ShortOpenWidth");
                double soh = chunk.GetDouble("ShortOpenHeight");
                it.ShortOpenWidth = sow > 0 ? sow : it.ShortOpenWidth;
                it.ShortOpenHeight = soh > 0 ? soh : it.ShortOpenHeight;

                double th = chunk.GetDouble("Thickness");
                it.Thickness = th > 0 ? th : it.Thickness;

                it.OffsetZ = chunk.GetDouble("OffsetZ");
                it.OffsetXY = chunk.GetDouble("OffsetXY");

                _instances.Add(it);
            }

            return base.Read(reader);
        }

        protected override Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.PullBox.png"))
                {
                    return stream != null ? new Bitmap(stream) : null;
                }
            }
        }

        public override Guid ComponentGuid => new Guid("8A3C6B7E-CC3B-4AC8-9C34-5A0D1D2B8E11");
    }
}