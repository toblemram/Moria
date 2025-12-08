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
    /// Places multiple circular pipes inside the under-shelf box and
    /// makes them follow the box along the tunnel (including bays).
    /// </summary>
    public class GH_TunnelBoxPipes : GH_Component
    {
        public GH_TunnelBoxPipes()
            : base(
                "Tunnel Box Pipes",
                "BoxPipes",
                "Creates multiple pipes inside the under-shelf box and makes them follow the box along the tunnel (including bays).",
                "Tunnel",
                "Details")
        { }

        // --------------------------------------------------------------------
        // INPUT
        // --------------------------------------------------------------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter(
                "Path",
                "P",
                "Tunnel centerline/path (same as used for the shelf/barrier component).",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "UnderShelfBox",
                "UB",
                "Box under the shelf (output from the barrier/shelf component).",
                GH_ParamAccess.item);

            p.AddNumberParameter(
                "PipeOffsetU",
                "U",
                "Horizontal offsets inside the box [m].\n" +
                "Measured from the local box corner with minimum X in each cross-section.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "PipeOffsetV",
                "V",
                "Vertical offsets inside the box [m].\n" +
                "Measured from the local box bottom (minimum Y in each cross-section).",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "PipeDiameters",
                "D",
                "Pipe diameters [m]. One value per pipe.",
                GH_ParamAccess.list);

            p.AddNumberParameter(
                "SectionStep",
                "ds",
                "Step along path for box sampling [m]. Smaller = smoother pipes but slower.",
                GH_ParamAccess.item,
                2.0);
        }

        // --------------------------------------------------------------------
        // OUTPUT
        // --------------------------------------------------------------------

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter(
                "Pipes",
                "P",
                "Pipe solids following the box along the tunnel.",
                GH_ParamAccess.list);

            p.AddPointParameter(
                "PipeCentersStart",
                "C0",
                "Pipe centres in the start cross-section.",
                GH_ParamAccess.list);

            p.AddTextParameter(
                "Info",
                "i",
                "Diagnostics/messages.",
                GH_ParamAccess.list);

            p.AddGeometryParameter(
                "Debug",
                "D",
                "Debug geometry (sampled centres, sections, etc.).",
                GH_ParamAccess.list);
        }

        // --------------------------------------------------------------------
        // MAIN SOLVE
        // --------------------------------------------------------------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            Brep box = null;
            var offsU = new List<double>();
            var offsV = new List<double>();
            var diams = new List<double>();
            double sectionStep = 2.0;

            da.GetData(0, ref path);
            da.GetData(1, ref box);
            da.GetDataList(2, offsU);
            da.GetDataList(3, offsV);
            da.GetDataList(4, diams);
            da.GetData(5, ref sectionStep);

            var info = new List<string>();
            var debug = new List<GeometryBase>();
            var pipes = new List<Brep>();
            var centersStart = new List<Point3d>();

            if (path == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Path is null.");
                return;
            }
            if (box == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "UnderShelfBox is null.");
                return;
            }

            int n = Math.Min(diams.Count, Math.Min(offsU.Count, offsV.Count));
            if (n == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No pipes defined (check U, V and D lists).");
                da.SetDataList(0, pipes);
                da.SetDataList(1, centersStart);
                da.SetDataList(2, info);
                da.SetDataList(3, debug);
                return;
            }

            if (sectionStep <= 0.1) sectionStep = 0.1;
            double tol = RhinoDoc.ActiveDoc != null
                ? RhinoDoc.ActiveDoc.ModelAbsoluteTolerance
                : 1e-3;

            double pathLen = path.GetLength();

            // ----------------------------------------------------------------
            // 1) Start-section, bare for å rapportere start-senterne
            // ----------------------------------------------------------------
            if (!path.LengthParameter(0.0, out double t0))
                t0 = path.Domain.Min;

            if (!path.PerpendicularFrameAt(t0, out Plane frame0))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Failed to get perpendicular frame at start of path.");
                return;
            }

            Curve[] cs0;
            Point3d[] pts0;
            if (!Intersection.BrepPlane(box, frame0, tol, out cs0, out pts0) ||
                cs0 == null || cs0.Length == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Can not intersect UnderShelfBox with start frame.");
                return;
            }

            var to2D_0 = Transform.PlaneToPlane(frame0, Plane.WorldXY);
            BoundingBox bb0_2D = BoundingBox.Empty;
            foreach (var c in cs0)
            {
                Curve c2 = c.DuplicateCurve();
                c2.Transform(to2D_0);
                var bb = c2.GetBoundingBox(true);
                if (!bb0_2D.IsValid) bb0_2D = bb;
                else bb0_2D.Union(bb);
            }

            double baseX0 = bb0_2D.Min.X;
            double baseY0 = bb0_2D.Min.Y;

            for (int i = 0; i < n; i++)
            {
                double u = offsU[i];
                double v = offsV[i];
                double d = diams[i];
                if (d <= 0.0) continue;

                double cx = baseX0 + u;
                double cy = baseY0 + v;

                Point3d p2D = new Point3d(cx, cy, 0.0);
                p2D.Transform(Transform.PlaneToPlane(Plane.WorldXY, frame0));
                centersStart.Add(p2D);
                debug.Add(new Rhino.Geometry.Point(p2D));
            }

            // ----------------------------------------------------------------
            // 2) For hvert rør: sample boksen langs hele path og bygg pipe
            // ----------------------------------------------------------------
            for (int i = 0; i < n; i++)
            {
                double u = offsU[i];
                double v = offsV[i];
                double d = diams[i];
                double r = 0.5 * d;

                if (d <= 0.0)
                {
                    info.Add($"Pipe {i}: diameter <= 0, skipped.");
                    continue;
                }

                var centerPts = new List<Point3d>();

                // Sampler langs tunnelen
                for (double s = 0.0; s <= pathLen + 0.5 * sectionStep; s += sectionStep)
                {
                    double clampedS = Math.Min(s, pathLen);

                    if (!path.LengthParameter(clampedS, out double tS))
                        continue;

                    if (!path.PerpendicularFrameAt(tS, out Plane frameS))
                        continue;

                    Curve[] cs;
                    Point3d[] dummy;
                    if (!Intersection.BrepPlane(box, frameS, tol, out cs, out dummy) ||
                        cs == null || cs.Length == 0)
                    {
                        continue;
                    }

                    // 2D-seksjon for denne posisjonen
                    var to2D = Transform.PlaneToPlane(frameS, Plane.WorldXY);
                    BoundingBox bb2D = BoundingBox.Empty;
                    foreach (var c in cs)
                    {
                        Curve c2 = c.DuplicateCurve();
                        c2.Transform(to2D);
                        var bb = c2.GetBoundingBox(true);
                        if (!bb2D.IsValid) bb2D = bb;
                        else bb2D.Union(bb);
                    }

                    if (!bb2D.IsValid)
                        continue;

                    double baseX = bb2D.Min.X;
                    double baseY = bb2D.Min.Y;

                    double cx = baseX + u;
                    double cy = baseY + v;

                    Point3d p2D = new Point3d(cx, cy, 0.0);
                    p2D.Transform(Transform.PlaneToPlane(Plane.WorldXY, frameS));

                    centerPts.Add(p2D);
                    debug.Add(new Rhino.Geometry.Point(p2D));
                }

                if (centerPts.Count < 2)
                {
                    info.Add($"Pipe {i}: too few sample points to build a centre-line.");
                    continue;
                }

                // Lag en glatt midtlinje
                Curve centerCurve = Curve.CreateInterpolatedCurve(centerPts, 3);
                if (centerCurve == null)
                {
                    info.Add($"Pipe {i}: failed to build interpolated centre curve.");
                    continue;
                }

                // Start- og endtverrsnitt for selve røret
                centerCurve.Domain = new Interval(0.0, 1.0);
                centerCurve.PerpendicularFrameAt(0.0, out Plane sPlane);
                centerCurve.PerpendicularFrameAt(1.0, out Plane ePlane);

                var circleStart = new Circle(sPlane, centerCurve.PointAtStart, r);
                var circleEnd = new Circle(ePlane, centerCurve.PointAtEnd, r);

                Curve cStart = circleStart.ToNurbsCurve();
                Curve cEnd = circleEnd.ToNurbsCurve();

                var sweep = new SweepOneRail
                {
                    SweepTolerance = tol,
                    AngleToleranceRadians = RhinoMath.ToRadians(1.0)
                };

                Brep[] pipeBreps = sweep.PerformSweep(centerCurve, new Curve[] { cStart, cEnd });
                if (pipeBreps == null || pipeBreps.Length == 0)
                {
                    info.Add($"Pipe {i}: sweep failed.");
                    continue;
                }

                Brep pipe = pipeBreps[0];
                if (!pipe.IsSolid)
                {
                    Brep capped = pipe.CapPlanarHoles(tol);
                    if (capped != null) pipe = capped;
                }

                pipes.Add(pipe);
            }

            // ----------------------------------------------------------------
            // OUTPUT
            // ----------------------------------------------------------------

            da.SetDataList(0, pipes);
            da.SetDataList(1, centersStart);
            da.SetDataList(2, info);
            da.SetDataList(3, debug);
        }

        // --------------------------------------------------------------------
        // MISC
        // --------------------------------------------------------------------

        protected override Bitmap Icon => null;

        public override Guid ComponentGuid => new Guid("7F6A7E0B-2F0D-4B0A-9B39-9B9F0A0CF9A2");
    }
}
