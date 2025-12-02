using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Generates concrete wall panels along tunnel walls based on a computed
    /// outer contact polyline per side (left/right), including emergency bay
    /// widening. 
    /// 
    /// Steps:
    ///   1. Sample the tunnel path and compute interpolated 2D profiles 
    ///      (base/bay) to obtain the wall "contact line" for each side.
    ///   2. Build a 2D curved plate profile based on the constant wall radius Rv
    ///      and center shift Yv from the T-profile parameters.
    ///   3. For each plate:
    ///         - place the 2D profile perpendicular to the wall polyline
    ///         - extrude the profile along the tangent of the wall polyline 
    ///           for a fixed panel length L
    ///   4. Spacing is applied along the wall polyline, independently per side.
    /// </summary>
    public class GH_TunnelConcretePanels : GH_Component
    {
        public GH_TunnelConcretePanels()
            : base("Tunnel.ConcretePanels", "Panels",
                   "Generates curved concrete wall panels following the tunnel wall, including emergency bay widenings.",
                   "Tunnel", "Details")
        { }

        // ---------------- INPUTS ----------------
        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            // Base profile FIRST
            p.AddTextParameter("BaseProfile", "Tbase",
                "Base T-profile: T5.5–T14.", GH_ParamAccess.item, "T10.5");

            // Path SECOND
            p.AddCurveParameter("Path", "P",
                "Tunnel centerline/path.", GH_ParamAccess.item);

            p.AddNumberParameter("BayStationsRight", "SR",
                "Stations (m) for emergency bays on the RIGHT side.", GH_ParamAccess.list);

            p.AddNumberParameter("BayStationsLeft", "SL",
                "Stations (m) for emergency bays on the LEFT side.", GH_ParamAccess.list);

            p.AddNumberParameter("Panel Length L", "L",
                "Panel length along wall [m].", GH_ParamAccess.item, 2.0);

            p.AddNumberParameter("Panel Height H", "H",
                "Panel height measured along arc [m].", GH_ParamAccess.item, 2.4);

            p.AddNumberParameter("Panel Thickness t", "t",
                "Wall panel thickness inward [m].", GH_ParamAccess.item, 0.25);

            p.AddNumberParameter("Spacing", "S",
                "Gap between panels along wall [m].", GH_ParamAccess.item, 0.0);

            p.AddBooleanParameter("MakeLeft", "Left",
                "Generate panels on left side.", GH_ParamAccess.item, true);

            p.AddBooleanParameter("MakeRight", "Right",
                "Generate panels on right side.", GH_ParamAccess.item, true);
        }

        // ---------------- OUTPUTS ----------------
        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Panels", "B",
                "Concrete wall panel Breps.", GH_ParamAccess.list);

            p.AddCurveParameter("LeftWallCurve", "CL",
                "Computed wall contact polyline on LEFT side.", GH_ParamAccess.item);

            p.AddCurveParameter("RightWallCurve", "CR",
                "Computed wall contact polyline on RIGHT side.", GH_ParamAccess.item);

            p.AddTextParameter("Info", "i",
                "Debug / status information.", GH_ParamAccess.list);
        }

        // ---------------- INTERNAL TYPES ----------------

        private class BayDef
        {
            public double SCenter;
            public int Side; // +1 = right, -1 = left
            public double S0, S1, S2, S3;
        }

        private class SegmentDef
        {
            public double SStart, SEnd;
            public Curve ShapeStart2D, ShapeEnd2D;
        }

        // ---------------- MAIN SOLVE ----------------
        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string baseProfile = "T10.5";
            var baysRight = new List<double>();
            var baysLeft = new List<double>();
            double L = 2.0, H = 2.4, t = 0.25, spacing = 0.0;
            bool makeLeft = true, makeRight = true;

            da.GetData(0, ref baseProfile);
            if (!da.GetData(1, ref path) || path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Missing tunnel path.");
                return;
            }

            da.GetDataList(2, baysRight);
            da.GetDataList(3, baysLeft);
            da.GetData(4, ref L);
            da.GetData(5, ref H);
            da.GetData(6, ref t);
            da.GetData(7, ref spacing);
            da.GetData(8, ref makeLeft);
            da.GetData(9, ref makeRight);

            var info = new List<string>();
            var panels = new List<Brep>();

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double pathLen = path.GetLength();
            if (pathLen <= tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path too short.");
                return;
            }

            baseProfile = baseProfile.Replace(",", ".").ToUpperInvariant();

            // --- Load T-profile parameters
            if (!ProfileType.Profiles.TryGetValue(baseProfile, out var basePar))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Invalid base profile '{baseProfile}'.");
                return;
            }

            string bayProfile = ResolveBayProfileType(baseProfile);
            if (!ProfileType.Profiles.TryGetValue(bayProfile, out var bayPar))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Bay profile '{bayProfile}' not found.");
                return;
            }

            info.Add($"Base profile: {baseProfile}");
            info.Add($"Bay profile  : {bayProfile}");

            // --- Build base and bay T-profiles in WorldXY
            if (!BuildProfileWorld(baseProfile, basePar, tol, out Curve base2D, out string err1))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, err1);
                return;
            }

            if (!BuildProfileWorld(bayProfile, bayPar, tol, out Curve bay2D, out string err2))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, err2);
                return;
            }

            // Horizontal shift logic matching other components
            BoundingBox bbBase = base2D.GetBoundingBox(true);
            BoundingBox bbBay = bay2D.GetBoundingBox(true);
            double xMinBase = bbBase.Min.X, xMaxBase = bbBase.Max.X;
            double xMinBay = bbBay.Min.X, xMaxBay = bbBay.Max.X;

            double dxRight = xMinBase - xMinBay;
            double dxLeft = xMaxBase - xMaxBay;

            Curve bayRight2D = bay2D.DuplicateCurve();
            bayRight2D.Transform(Transform.Translation(dxRight, 0, 0));

            Curve bayLeft2D = bay2D.DuplicateCurve();
            bayLeft2D.Transform(Transform.Translation(dxLeft, 0, 0));

            // --- Build bay definitions and segments
            var bayDefs = BuildBayDefs(baysRight, baysLeft, pathLen, info);
            var segments = BuildSegments(bayDefs, base2D, bayRight2D, bayLeft2D, pathLen);

            // --- Build left/right wall polylines
            PolylineCurve leftWall = null;
            PolylineCurve rightWall = null;

            int samples = 200;
            if (makeLeft)
                leftWall = BuildWallPolyline(path, segments, isLeft: true, samples, tol);
            if (makeRight)
                rightWall = BuildWallPolyline(path, segments, isLeft: false, samples, tol);

            // --- Build curved 2D plate profile (arc sector)
            Curve plate2D = BuildPlateSection2D(H, t, basePar.Rv, basePar.Yv, tol);
            if (plate2D == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Failed to build 2D plate profile.");
                return;
            }

            double pitch = L + spacing;

            // --- Generate panels on each side
            if (makeLeft && leftWall != null)
                GeneratePanelsOnSide(leftWall, true, L, pitch, plate2D, tol, panels);

            if (makeRight && rightWall != null)
                GeneratePanelsOnSide(rightWall, false, L, pitch, plate2D, tol, panels);

            info.Add($"Total panels: {panels.Count}");

            // --- Output
            da.SetDataList(0, panels);
            da.SetData(1, leftWall);
            da.SetData(2, rightWall);
            da.SetDataList(3, info);
        }

        // ---------------- PANEL GENERATION ----------------
        private void GeneratePanelsOnSide(
            Curve wallCurve,
            bool isLeft,
            double L,
            double pitch,
            Curve plate2D,
            double tol,
            List<Brep> panels)
        {
            double wallLen = wallCurve.GetLength();
            if (wallLen <= tol) return;

            int n = (int)Math.Floor((wallLen - L) / pitch + 1e-9);
            if (n < 0) n = 0;

            for (int i = 0; i <= n; i++)
            {
                double s0 = i * pitch;
                double s1 = s0 + L;
                if (s1 > wallLen + 1e-6) break;

                if (!wallCurve.LengthParameter(s0, out double t0)) continue;
                if (!wallCurve.LengthParameter(s1, out double t1)) continue;

                Point3d B0 = wallCurve.PointAt(t0);
                Point3d B1 = wallCurve.PointAt(t1);

                Vector3d forward = B1 - B0;
                if (!forward.Unitize() || forward.IsTiny()) continue;

                // Global vertical
                Vector3d up = Vector3d.ZAxis;
                if (!up.Unitize()) continue;

                // Side direction in XY-plane
                Vector3d side = Vector3d.CrossProduct(up, forward);
                if (!side.Unitize()) continue;

                // Inward direction
                Vector3d inward = isLeft ? -side : side;
                inward.Unitize();

                // Local coordinate plane for the 2D profile:
                //   origin = inner bottom corner at B0
                //   X axis = inward
                //   Y axis = up (panel curvature direction)
                Plane platePlane = new Plane(B0, inward, up);

                // Position 2D profile
                Curve sec3d = plate2D.DuplicateCurve();
                sec3d.Transform(Transform.PlaneToPlane(Plane.WorldXY, platePlane));

                // Extrude along wall tangent
                Vector3d extrude = forward * L;
                Surface srf = Surface.CreateExtrusion(sec3d, extrude);
                if (srf == null) continue;

                Brep b = srf.ToBrep().CapPlanarHoles(tol);
                if (b != null) panels.Add(b);
            }
        }

        // ---------------- BUILD WALL POLYLINE ----------------
        private PolylineCurve BuildWallPolyline(
            Curve path,
            List<SegmentDef> segments,
            bool isLeft,
            int samples,
            double tol)
        {
            var pts = new List<Point3d>();
            double pathLen = path.GetLength();

            for (int i = 0; i <= samples; i++)
            {
                double sPath = pathLen * i / (double)samples;

                if (!path.LengthParameter(sPath, out double tPath)) continue;
                if (!path.PerpendicularFrameAt(tPath, out Plane frame)) continue;

                if (!FindSegment(segments, sPath, out SegmentDef seg, out double uSeg)) continue;

                Curve prof2D = InterpolateSection(seg.ShapeStart2D, seg.ShapeEnd2D, uSeg, tol);
                if (!(prof2D is PolyCurve pc) || pc.SegmentCount < 4) continue;

                // 0 = left wall segment, 2 = right wall segment
                Curve wall2D = isLeft ? pc.SegmentCurve(0) : pc.SegmentCurve(2);

                Point3d p0 = wall2D.PointAtStart;
                Point3d p1 = wall2D.PointAtEnd;
                Point3d bottom = (p0.Y <= p1.Y) ? p0 : p1;

                Point3d worldPt = bottom;
                worldPt.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame));

                pts.Add(worldPt);
            }

            if (pts.Count < 2) return null;

            return new PolylineCurve(pts);
        }

        // ---------------- BUILD CURVED PANEL SECTION ----------------
        /// <summary>
        /// Builds a curved 2D panel cross-section using:
        ///   - Circle radius Rv
        ///   - Circle center (0, Yv)
        ///   - Starting at arc position where circle intersects y = 0
        ///   - Height H measured along the arc
        ///   - Thickness t = radial expansion
        ///
        /// After generating both arcs, the inner bottom is moved to (0,0)
        /// so the 2D profile can be mapped to the wall polyline naturally.
        /// </summary>
        private Curve BuildPlateSection2D(double H, double t, double Rv, double Yv, double tol)
        {
            if (H <= tol || t <= tol)
                return null;

            // Solve y = Yv + Rv*sin(theta) = 0  →  sin(theta0) = -Yv/Rv
            double ratio = -Yv / Rv;
            if (ratio < -1.0) ratio = -1.0;
            if (ratio > 1.0) ratio = 1.0;

            double theta0 = Math.Asin(ratio);
            double theta1 = theta0 + (H / Rv);

            int n = 32;
            var innerPts = new List<Point3d>();
            var outerPts = new List<Point3d>();

            for (int i = 0; i <= n; i++)
            {
                double f = (double)i / n;
                double th = theta0 + (theta1 - theta0) * f;

                double cos = Math.Cos(th);
                double sin = Math.Sin(th);

                // inner arc
                double ix = Rv * cos;
                double iy = Yv + Rv * sin;

                // outer arc
                double Rout = Rv + t;
                double ox = Rout * cos;
                double oy = Yv + Rout * sin;

                innerPts.Add(new Point3d(ix, iy, 0.0));
                outerPts.Add(new Point3d(ox, oy, 0.0));
            }

            // Move inner bottom to (0,0)
            Vector3d shift = -(Vector3d)innerPts[0];

            for (int i = 0; i < innerPts.Count; i++)
                innerPts[i] = innerPts[i] + shift;
            for (int i = 0; i < outerPts.Count; i++)
                outerPts[i] = outerPts[i] + shift;

            // Build closed polyline: inner bottom→top, outer top→bottom
            outerPts.Reverse();

            var ptsAll = new List<Point3d>();
            ptsAll.AddRange(innerPts);
            ptsAll.AddRange(outerPts);
            ptsAll.Add(innerPts[0]); // close

            return new Polyline(ptsAll).ToNurbsCurve();
        }

        // ---------------- BAY + SEGMENT SETUP ----------------
        private List<BayDef> BuildBayDefs(
            List<double> baysRight,
            List<double> baysLeft,
            double pathLength,
            List<string> info)
        {
            const double ramp = 30.0;
            const double flat = 30.0;

            var list = new List<BayDef>();

            foreach (double s in baysRight)
            {
                list.Add(new BayDef
                {
                    SCenter = s,
                    Side = +1,
                    S0 = s - (ramp + flat * 0.5),
                    S1 = s - (flat * 0.5),
                    S2 = s + (flat * 0.5),
                    S3 = s + (ramp + flat * 0.5)
                });
            }

            foreach (double s in baysLeft)
            {
                list.Add(new BayDef
                {
                    SCenter = s,
                    Side = -1,
                    S0 = s - (ramp + flat * 0.5),
                    S1 = s - (flat * 0.5),
                    S2 = s + (flat * 0.5),
                    S3 = s + (ramp + flat * 0.5)
                });
            }

            list.Sort((a, b) => a.SCenter.CompareTo(b.SCenter));

            foreach (var b in list)
            {
                if (b.S0 < 0)
                {
                    info?.Add($"Bay at s={b.SCenter:0.###}: S0 < 0 → set to 0.");
                    b.S0 = 0;
                }
                if (b.S3 > pathLength)
                {
                    info?.Add($"Bay at s={b.SCenter:0.###}: S3 > path length → trimmed.");
                    b.S3 = pathLength;
                }
            }

            return list;
        }

        private List<SegmentDef> BuildSegments(
            List<BayDef> bays,
            Curve base2D,
            Curve bayRight2D,
            Curve bayLeft2D,
            double pathLength)
        {
            var segments = new List<SegmentDef>();
            double cur = 0.0;
            const double eps = 1e-6;

            foreach (var b in bays)
            {
                if (b.S0 > cur + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = cur,
                        SEnd = b.S0,
                        ShapeStart2D = base2D,
                        ShapeEnd2D = base2D
                    });
                }

                Curve bay2D = (b.Side > 0) ? bayRight2D : bayLeft2D;

                if (b.S1 > b.S0 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S0,
                        SEnd = b.S1,
                        ShapeStart2D = base2D,
                        ShapeEnd2D = bay2D
                    });
                }

                if (b.S2 > b.S1 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S1,
                        SEnd = b.S2,
                        ShapeStart2D = bay2D,
                        ShapeEnd2D = bay2D
                    });
                }

                if (b.S3 > b.S2 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S2,
                        SEnd = b.S3,
                        ShapeStart2D = bay2D,
                        ShapeEnd2D = base2D
                    });
                }

                cur = Math.Max(cur, b.S3);
            }

            if (cur < pathLength - eps)
            {
                segments.Add(new SegmentDef
                {
                    SStart = cur,
                    SEnd = pathLength,
                    ShapeStart2D = base2D,
                    ShapeEnd2D = base2D
                });
            }

            if (bays.Count == 0 && segments.Count == 0)
            {
                segments.Add(new SegmentDef
                {
                    SStart = 0.0,
                    SEnd = pathLength,
                    ShapeStart2D = base2D,
                    ShapeEnd2D = base2D
                });
            }

            return segments;
        }

        private bool FindSegment(List<SegmentDef> segments, double sPath,
                                 out SegmentDef seg, out double u)
        {
            seg = null;
            u = 0.0;
            const double eps = 1e-6;

            foreach (var sg in segments)
            {
                if (sPath + eps < sg.SStart || sPath - eps > sg.SEnd)
                    continue;

                seg = sg;
                double len = sg.SEnd - sg.SStart;
                if (len <= eps) return true;

                u = (sPath - sg.SStart) / len;
                u = Math.Max(0.0, Math.Min(1.0, u));
                return true;
            }

            return false;
        }

        // ---------------- PROFILE GENERATION + INTERPOLATION ----------------
        private bool BuildProfileWorld(string type, ProfileType.ProfileParameters par,
                                       double tol, out Curve profile, out string error)
        {
            profile = null;
            error = null;

            PolyCurve poly;
            var dummySeg = new List<Curve>();
            var dummyDbg = new List<GeometryBase>();
            bool ok;

            if (ProfileType.IsLowRoof(type))
            {
                ok = LowRoofProfileBuilder.Build(type, par, true, tol,
                                                 out poly, out dummySeg, out dummyDbg, out error);
            }
            else
            {
                ok = StandardProfileBuilder.Build(par, true, tol,
                                                  out poly, out dummySeg, out dummyDbg, out error);
            }

            if (!ok || poly == null)
                return false;

            // translate bottom midpoint to origin
            Curve bottom = poly.SegmentCurve(poly.SegmentCount - 1);
            Point3d b0 = bottom.PointAtStart;
            Point3d b1 = bottom.PointAtEnd;
            Point3d mid = 0.5 * (b0 + b1);

            Transform toOrigin = Transform.Translation(-mid.X, -mid.Y, -mid.Z);
            poly.Transform(toOrigin);

            profile = poly;
            return true;
        }

        private string ResolveBayProfileType(string baseProfile)
        {
            switch (baseProfile)
            {
                case "T5.5": return "T8.5";
                case "T7.5": return "T10.5";
                case "T8.5": return "T12.5";
                default: return "T14";
            }
        }

        private Curve InterpolateSection(Curve start, Curve end, double u, double tol)
        {
            u = Math.Max(0.0, Math.Min(1.0, u));

            if (u <= tol) return start.DuplicateCurve();
            if (u >= 1.0 - tol) return end.DuplicateCurve();
            if (!(start is PolyCurve ps) || !(end is PolyCurve pe))
                return null;
            if (ps.SegmentCount != pe.SegmentCount)
                return null;

            var res = new PolyCurve();
            int segCount = ps.SegmentCount;
            int samplesPerSeg = 24;

            for (int i = 0; i < segCount; i++)
            {
                Curve a = ps.SegmentCurve(i);
                Curve b = pe.SegmentCurve(i);

                var pts = new List<Point3d>();
                for (int k = 0; k <= samplesPerSeg; k++)
                {
                    double tn = (double)k / samplesPerSeg;
                    Point3d pa = a.PointAtNormalizedLength(tn);
                    Point3d pb = b.PointAtNormalizedLength(tn);
                    pts.Add((1.0 - u) * pa + u * pb);
                }

                res.AppendSegment(new Polyline(pts).ToNurbsCurve());
            }

            res.MakeClosed(tol);
            return res;
        }

        // ---------------- METADATA ----------------
        public override Guid ComponentGuid => new Guid("A3A5F7DA-40E6-4EBE-A8AE-68336D2D2F31");
        protected override System.Drawing.Bitmap Icon => null;
    }
}
