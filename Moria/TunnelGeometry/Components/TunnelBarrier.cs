using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Moria.TunnelGeometry.Components
{
    public class GH_TunnelBarrier : GH_Component
    {
        public GH_TunnelBarrier()
            : base(
                "Tunnel Barrier + Shelf",
                "BarrierShelf",
                "Creates a concrete barrier and upper shelf along one tunnel wall, following T-profiles and emergency bays.",
                "Tunnel",
                "Details")
        { }

        // --------------------------------------------------------------------
        // INPUT
        // --------------------------------------------------------------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter("Path", "P", "Tunnel centerline/path", GH_ParamAccess.item);

            p.AddTextParameter(
                "BaseProfileType",
                "T",
                "Base tunnel T-profile, e.g. T10.5, T12.5, T14.",
                GH_ParamAccess.item,
                "T10.5");

            p.AddNumberParameter(
                "BayStationsRight",
                "SR",
                "Bay stations on right side.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "BayStationsLeft",
                "SL",
                "Bay stations on left side.",
                GH_ParamAccess.list);

            p.AddBooleanParameter(
                "RightSide",
                "R",
                "True = place barrier on right side. False = left side.",
                GH_ParamAccess.item,
                true);

            p.AddNumberParameter(
                "BarrierOffset",
                "BO",
                "Vertical offset of barrier base from tunnel bottom [m]",
                GH_ParamAccess.item,
                0.30);

            p.AddNumberParameter(
                "BarrierHeight",
                "BH",
                "Barrier height [m]",
                GH_ParamAccess.item,
                1.00);

            p.AddNumberParameter(
                "BarrierBaseThickness",
                "BB",
                "Barrier thickness at bottom [m]",
                GH_ParamAccess.item,
                0.30);

            p.AddNumberParameter(
                "BarrierTopThickness",
                "BT",
                "Barrier thickness at top [m]",
                GH_ParamAccess.item,
                0.20);

            p.AddNumberParameter(
                "GapBetweenShelfAndBarrier",
                "GW",
                "Gap between front of shelf and back of barrier [m]",
                GH_ParamAccess.item,
                0.05);

            p.AddNumberParameter(
                "ShelfWidth",
                "SW",
                "Shelf width measured from wall toward roadway [m]",
                GH_ParamAccess.item,
                0.40);

            p.AddNumberParameter(
                "ShelfThickness",
                "ST",
                "Vertical thickness of shelf [m]",
                GH_ParamAccess.item,
                0.15);
        }

        // --------------------------------------------------------------------
        // OUTPUT
        // --------------------------------------------------------------------

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Barrier", "B", "Barrier geometry", GH_ParamAccess.item);
            p.AddBrepParameter("Shelf", "S", "Shelf geometry", GH_ParamAccess.item);
            p.AddCurveParameter("EdgeLine", "E", "Edge line (top of shelf)", GH_ParamAccess.item);
            p.AddTextParameter("Info", "i", "Diagnostics", GH_ParamAccess.list);
            p.AddGeometryParameter("Debug", "D", "Debug geometry", GH_ParamAccess.list);
        }

        // --------------------------------------------------------------------
        // INTERNAL TYPES
        // --------------------------------------------------------------------

        private class BayDef
        {
            public double SCenter;
            public int Side;   // +1 right, -1 left
            public double S0;
            public double S1;
            public double S2;
            public double S3;
        }

        private class SegmentDef
        {
            public double SStart;
            public double SEnd;

            public Curve BarrierStart2D;
            public Curve BarrierEnd2D;

            public Curve ShelfStart2D;
            public Curve ShelfEnd2D;

            public Point3d EdgeStart2D;
            public Point3d EdgeEnd2D;
        }

        // --------------------------------------------------------------------
        // MAIN SOLVE
        // --------------------------------------------------------------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string baseProfileType = "T10.5";
            var stationsRight = new List<double>();
            var stationsLeft = new List<double>();
            bool rightSide = true;

            double barrierOffset = 0.30;
            double barrierHeight = 1.00;
            double barrierBaseThickness = 0.30;
            double barrierTopThickness = 0.20;
            double gapShelfBarrier = 0.05;
            double shelfWidth = 0.40;
            double shelfThickness = 0.15;

            da.GetData(0, ref path);
            da.GetData(1, ref baseProfileType);
            da.GetDataList(2, stationsRight);
            da.GetDataList(3, stationsLeft);
            da.GetData(4, ref rightSide);
            da.GetData(5, ref barrierOffset);
            da.GetData(6, ref barrierHeight);
            da.GetData(7, ref barrierBaseThickness);
            da.GetData(8, ref barrierTopThickness);
            da.GetData(9, ref gapShelfBarrier);
            da.GetData(10, ref shelfWidth);
            da.GetData(11, ref shelfThickness);

            var info = new List<string>();
            var debug = new List<GeometryBase>();

            if (path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path is null");
                return;
            }

            double tol = RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
            double pathLength = path.GetLength();

            // ------------------------
            // BAY DEFINITIONS
            // ------------------------

            const double rampLen = 30.0;
            const double flatLen = 30.0;

            var bays = new List<BayDef>();

            if (rightSide)
            {
                foreach (var s in stationsRight)
                {
                    bays.Add(new BayDef
                    {
                        SCenter = s,
                        Side = +1,
                        S0 = s - (rampLen + flatLen * 0.5),
                        S1 = s - (flatLen * 0.5),
                        S2 = s + (flatLen * 0.5),
                        S3 = s + (rampLen + flatLen * 0.5)
                    });
                }
            }
            else
            {
                foreach (var s in stationsLeft)
                {
                    bays.Add(new BayDef
                    {
                        SCenter = s,
                        Side = -1,
                        S0 = s - (rampLen + flatLen * 0.5),
                        S1 = s - (flatLen * 0.5),
                        S2 = s + (flatLen * 0.5),
                        S3 = s + (rampLen + flatLen * 0.5)
                    });
                }
            }

            bays.Sort((a, b) => a.SCenter.CompareTo(b.SCenter));

            foreach (var b in bays)
            {
                if (b.S0 < 0) b.S0 = 0;
                if (b.S3 > pathLength) b.S3 = pathLength;
            }

            // ------------------------
            // BUILD BASE + BAY PROFILES
            // ------------------------

            baseProfileType = baseProfileType.ToUpper().Replace(",", ".");
            if (!ProfileType.Profiles.TryGetValue(baseProfileType, out var baseParams))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Invalid base profile: {baseProfileType}");
                return;
            }

            string bayType = ResolveBayProfileType(baseProfileType);
            if (!ProfileType.Profiles.TryGetValue(bayType, out var bayParams))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Invalid bay profile: {bayType}");
                return;
            }

            BuildProfileWorld(baseProfileType, baseParams, tol, out Curve baseProf2D, out _);
            BuildProfileWorld(bayType, bayParams, tol, out Curve bayProf2D, out _);

            // Shift bay profile to attach correctly to base profile (same as EmergencyBay)
            var bbBase = baseProf2D.GetBoundingBox(true);
            var bbBay = bayProf2D.GetBoundingBox(true);

            double dxRight = bbBase.Min.X - bbBay.Min.X;
            double dxLeft = bbBase.Max.X - bbBay.Max.X;

            Curve bayRight = bayProf2D.DuplicateCurve();
            bayRight.Transform(Transform.Translation(dxRight, 0, 0));

            Curve bayLeft = bayProf2D.DuplicateCurve();
            bayLeft.Transform(Transform.Translation(dxLeft, 0, 0));

            Curve baySide2D = rightSide ? bayRight : bayLeft;

            debug.Add(baseProf2D.DuplicateCurve());
            debug.Add(baySide2D.DuplicateCurve());

            // ------------------------
            // BUILD BARRIER/SHELF 2D
            // ------------------------

            BuildBarrierAndShelf2D(
                baseProf2D,
                rightSide,
                barrierOffset,
                barrierHeight,
                barrierBaseThickness,
                barrierTopThickness,
                gapShelfBarrier,
                shelfWidth,
                shelfThickness,
                tol,
                out Curve barrierBase2D,
                out Curve shelfBase2D,
                out Point3d edgeBase2D,
                out string err1);

            BuildBarrierAndShelf2D(
                baySide2D,
                rightSide,
                barrierOffset,
                barrierHeight,
                barrierBaseThickness,
                barrierTopThickness,
                gapShelfBarrier,
                shelfWidth,
                shelfThickness,
                tol,
                out Curve barrierBay2D,
                out Curve shelfBay2D,
                out Point3d edgeBay2D,
                out string err2);

            debug.Add(barrierBase2D.DuplicateCurve());
            debug.Add(shelfBase2D.DuplicateCurve());
            debug.Add(barrierBay2D.DuplicateCurve());
            debug.Add(shelfBay2D.DuplicateCurve());
            debug.Add(new Rhino.Geometry.Point(edgeBase2D));
            debug.Add(new Rhino.Geometry.Point(edgeBay2D));

            // ------------------------
            // BUILD SWEEP SEGMENTS
            // ------------------------

            var segments = new List<SegmentDef>();
            double currentS = 0;
            double eps = 1e-6;

            foreach (var b in bays)
            {
                if (b.S0 > currentS + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = currentS,
                        SEnd = b.S0,

                        BarrierStart2D = barrierBase2D,
                        BarrierEnd2D = barrierBase2D,
                        ShelfStart2D = shelfBase2D,
                        ShelfEnd2D = shelfBase2D,
                        EdgeStart2D = edgeBase2D,
                        EdgeEnd2D = edgeBase2D
                    });
                }

                if (b.S1 > b.S0 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S0,
                        SEnd = b.S1,

                        BarrierStart2D = barrierBase2D,
                        BarrierEnd2D = barrierBay2D,
                        ShelfStart2D = shelfBase2D,
                        ShelfEnd2D = shelfBay2D,
                        EdgeStart2D = edgeBase2D,
                        EdgeEnd2D = edgeBay2D
                    });
                }

                if (b.S2 > b.S1 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S1,
                        SEnd = b.S2,

                        BarrierStart2D = barrierBay2D,
                        BarrierEnd2D = barrierBay2D,
                        ShelfStart2D = shelfBay2D,
                        ShelfEnd2D = shelfBay2D,
                        EdgeStart2D = edgeBay2D,
                        EdgeEnd2D = edgeBay2D
                    });
                }

                if (b.S3 > b.S2 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S2,
                        SEnd = b.S3,

                        BarrierStart2D = barrierBay2D,
                        BarrierEnd2D = barrierBase2D,
                        ShelfStart2D = shelfBay2D,
                        ShelfEnd2D = shelfBase2D,
                        EdgeStart2D = edgeBay2D,
                        EdgeEnd2D = edgeBase2D
                    });
                }

                currentS = b.S3;
            }

            if (currentS < pathLength - eps)
            {
                segments.Add(new SegmentDef
                {
                    SStart = currentS,
                    SEnd = pathLength,

                    BarrierStart2D = barrierBase2D,
                    BarrierEnd2D = barrierBase2D,
                    ShelfStart2D = shelfBase2D,
                    ShelfEnd2D = shelfBase2D,
                    EdgeStart2D = edgeBase2D,
                    EdgeEnd2D = edgeBase2D
                });
            }

            // ------------------------
            // SWEEP
            // ------------------------

            var sweep = new SweepOneRail()
            {
                SweepTolerance = tol,
                AngleToleranceRadians = RhinoMath.ToRadians(1)
            };

            var barrierPieces = new List<Brep>();
            var shelfPieces = new List<Brep>();
            int segId = 0;

            foreach (var seg in segments)
            {
                double segLen = seg.SEnd - seg.SStart;
                if (segLen <= tol)
                {
                    info.Add($"Segment {segId}: too short");
                    continue;
                }

                if (!path.LengthParameter(seg.SStart, out double t0) ||
                    !path.LengthParameter(seg.SEnd, out double t1))
                {
                    info.Add($"Segment {segId}: failed LengthParameter");
                    continue;
                }

                double tMin = Math.Min(t0, t1);
                double tMax = Math.Max(t0, t1);

                Curve rail = path.DuplicateCurve().Trim(tMin, tMax);
                if (rail == null)
                {
                    info.Add($"Segment {segId}: rail trim failed");
                    continue;
                }

                if (!path.PerpendicularFrameAt(tMin, out Plane frame0) ||
                    !path.PerpendicularFrameAt(tMax, out Plane frame1))
                {
                    info.Add($"Segment {segId}: frame failed");
                    continue;
                }

                // Barrier
                Curve bS = seg.BarrierStart2D.DuplicateCurve();
                Curve bE = seg.BarrierEnd2D.DuplicateCurve();

                bS.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame0));
                bE.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame1));

                var bPiece = sweep.PerformSweep(rail, new Curve[] { bS, bE });
                if (bPiece != null && bPiece.Length > 0)
                    barrierPieces.Add(bPiece[0]);

                // Shelf
                Curve sS = seg.ShelfStart2D.DuplicateCurve();
                Curve sE = seg.ShelfEnd2D.DuplicateCurve();

                sS.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame0));
                sE.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame1));

                var sPiece = sweep.PerformSweep(rail, new Curve[] { sS, sE });
                if (sPiece != null && sPiece.Length > 0)
                    shelfPieces.Add(sPiece[0]);

                segId++;
            }

            var barrierJoined = JoinBrepsList(barrierPieces, tol);
            var shelfJoined = JoinBrepsList(shelfPieces, tol);

            barrierJoined = PrepareSolid(barrierJoined, tol);
            shelfJoined = PrepareSolid(shelfJoined, tol);

            // ------------------------
            // EDGE LINE
            // ------------------------

            Curve edgeLine = BuildEdgeLine(path, segments, tol);

            // ------------------------
            // OUTPUT
            // ------------------------

            da.SetData(0, barrierJoined);
            da.SetData(1, shelfJoined);
            da.SetData(2, edgeLine);
            da.SetDataList(3, info);
            da.SetDataList(4, debug);
        }

        // --------------------------------------------------------------------
        // PROFILE BUILDER
        // --------------------------------------------------------------------

        private string ResolveBayProfileType(string t)
        {
            switch (t)
            {
                case "T5.5": return "T8.5";
                case "T7.5": return "T10.5";
                case "T8.5": return "T12.5";
                default: return "T14";
            }
        }

        private void BuildProfileWorld(
            string type,
            ProfileType.ProfileParameters par,
            double tol,
            out Curve profile,
            out string error)
        {
            profile = null;
            error = null;

            PolyCurve poly;
            var segs = new List<Curve>();
            var dbg = new List<GeometryBase>();

            bool ok =
                ProfileType.IsLowRoof(type)
                ? LowRoofProfileBuilder.Build(type, par, true, tol, out poly, out segs, out dbg, out error)
                : StandardProfileBuilder.Build(par, true, tol, out poly, out segs, out dbg, out error);

            if (!ok)
                return;

            Curve bottom = poly.SegmentCurve(poly.SegmentCount - 1);
            var mid = 0.5 * (bottom.PointAtStart + bottom.PointAtEnd);
            poly.Transform(Transform.Translation(-mid.X, -mid.Y, -mid.Z));

            profile = poly;
        }

        // --------------------------------------------------------------------
        // NEW BARRIER + SHELF GEOMETRY (UPDATED AS AGREED!)
        // --------------------------------------------------------------------

        private static bool BuildBarrierAndShelf2D(
            Curve profileWorld,
            bool rightSide,
            double barrierOffset,
            double barrierHeight,
            double barrierBaseThickness,
            double barrierTopThickness,
            double gapShelfBarrier,
            double shelfWidth,
            double shelfThickness,
            double tol,
            out Curve barrier2D,
            out Curve shelf2D,
            out Point3d edgePoint2D,
            out string error)
        {
            error = null;
            barrier2D = null;
            shelf2D = null;
            edgePoint2D = new Point3d();

            if (!(profileWorld is PolyCurve poly))
            {
                error = "Profile 2D is not a polycurve";
                return false;
            }

            Curve wall = poly.SegmentCurve(rightSide ? 2 : 0);

            var horz = new LineCurve(
                new Point3d(-1000, barrierOffset, 0),
                new Point3d(+1000, barrierOffset, 0));

            var ccx = Intersection.CurveCurve(wall, horz, tol, tol);
            if (ccx.Count == 0)
            {
                error = "Failed to intersect wall at barrierOffset";
                return false;
            }

            double xWall = ccx[0].PointA.X;
            foreach (var cr in ccx)
            {
                if (rightSide && cr.PointA.X > xWall) xWall = cr.PointA.X;
                if (!rightSide && cr.PointA.X < xWall) xWall = cr.PointA.X;
            }

            double sign = rightSide ? +1.0 : -1.0;

            double y0 = barrierOffset;
            double y1 = barrierOffset + barrierHeight;

            double yShelfTop = y1;
            double yShelfBottom = y1 - shelfThickness;

            // Shelf near wall
            double xShelfBack = xWall;
            double xShelfFront = xShelfBack - sign * shelfWidth;

            // Gap between shelf and barrier
            double xBack = xShelfFront - sign * gapShelfBarrier;

            // Barrier front
            double xFrontBase = xBack - sign * barrierBaseThickness;
            double xFrontTop = xBack - sign * barrierTopThickness;

            // Barrier poly
            barrier2D = new PolylineCurve(new[]
            {
                new Point3d(xFrontBase, y0, 0),
                new Point3d(xFrontTop,  y1, 0),
                new Point3d(xBack,      y1, 0),
                new Point3d(xBack,      y0, 0),
                new Point3d(xFrontBase, y0, 0)
            });

            // Shelf poly (flat, adjacent to barrier)
            shelf2D = new PolylineCurve(new[]
            {
                new Point3d(xShelfFront, yShelfBottom, 0),
                new Point3d(xShelfBack,  yShelfBottom, 0),
                new Point3d(xShelfBack,  yShelfTop,    0),
                new Point3d(xShelfFront, yShelfTop,    0),
                new Point3d(xShelfFront, yShelfBottom, 0)
            });

            edgePoint2D = new Point3d(xShelfFront, yShelfTop, 0);

            return true;
        }

        // --------------------------------------------------------------------
        // EDGE LINE
        // --------------------------------------------------------------------

        private Curve BuildEdgeLine(Curve path, List<SegmentDef> segments, double tol)
        {
            double pathLen = path.GetLength();
            double step = 2.0;

            var pts = new List<Point3d>();

            for (double s = 0; s <= pathLen; s += step)
            {
                GetSegmentForS(segments, s, out var seg, out double u);
                Point3d e2D = Lerp(seg.EdgeStart2D, seg.EdgeEnd2D, u);

                if (!path.LengthParameter(s, out double t)) continue;
                if (!path.PerpendicularFrameAt(t, out Plane f)) continue;

                var p3 = e2D;
                p3.Transform(Transform.PlaneToPlane(Plane.WorldXY, f));

                pts.Add(p3);
            }

            if (pts.Count < 2)
                return null;

            return new PolylineCurve(pts);
        }

        private void GetSegmentForS(List<SegmentDef> segs, double s, out SegmentDef seg, out double u)
        {
            seg = segs[0];
            foreach (var sg in segs)
            {
                if (s >= sg.SStart && s <= sg.SEnd)
                {
                    seg = sg;
                    break;
                }
            }

            double len = seg.SEnd - seg.SStart;
            u = (len <= 0) ? 0 : (s - seg.SStart) / len;
            if (u < 0) u = 0;
            if (u > 1) u = 1;
        }

        private Point3d Lerp(Point3d a, Point3d b, double t)
        {
            return new Point3d(
                a.X + (b.X - a.X) * t,
                a.Y + (b.Y - a.Y) * t,
                a.Z + (b.Z - a.Z) * t
            );
        }

        // --------------------------------------------------------------------
        // UTILS
        // --------------------------------------------------------------------

        private Brep JoinBrepsList(List<Brep> breps, double tol)
        {
            if (breps.Count == 0) return null;
            if (breps.Count == 1) return breps[0];

            var joined = Brep.JoinBreps(breps, tol);
            return joined != null && joined.Length > 0 ? joined[0] : breps[0];
        }

        private Brep PrepareSolid(Brep b, double tol)
        {
            if (b == null) return null;

            if (!b.IsSolid)
            {
                Brep capped = b.CapPlanarHoles(tol);
                if (capped != null) b = capped;

                b.Repair(tol);
            }

            return b;
        }

        protected override Bitmap Icon => null;

        public override Guid ComponentGuid => new Guid("C8E1B752-33A5-4FDF-A335-8CC9FB11EF56");
    }
}
