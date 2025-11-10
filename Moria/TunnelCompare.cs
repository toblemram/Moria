using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;           // For ToolStripDropDown i GH-meny
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Tunnel.GH
{
    public class GH_TunnelBlastVolume_AB : GH_Component
    {
        // --- Felt for sist beregnede seksjonsdata (til WPF-viewer) ---
        private List<double> _lastStations;
        private List<Point3d> _lastSectionPoints;
        private List<double> _lastBlastAreas;
        private List<double> _lastScanAreas;
        private List<double> _lastDesignAreas;

        // NEW: 2D-profiler per snitt (XY) til WPF-tegning
        private List<List<Polyline>> _lastScanProfiles2D;
        private List<List<Polyline>> _lastDesignProfiles2D;
        private List<List<Polyline>> _lastBlastProfiles2D;

        public GH_TunnelBlastVolume_AB()
          : base(
              "Tunnel Blast Volume (Sections)",
              "BlastVolSec",
              "Tilnærmet sprengvolum ved å ta tverrsnitt langs senterlinjen mellom A–B og loft’e differansekurver (Design minus Scan).",
              "Tunnel",
              "Analysis")
        { }

        public override Guid ComponentGuid => new Guid("2E2C3D7E-5C8E-4A19-9D3F-0E6F1C8E7C11");
        protected override Bitmap Icon => null;
        public override GH_Exposure Exposure => GH_Exposure.secondary;

        // ----------------- INPUT / OUTPUT -----------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddMeshParameter("ScanMesh", "Scan",
                "Mesh av eksisterende tunnel (innvendig skann, åpen/lukket).",
                GH_ParamAccess.item);

            p.AddBrepParameter("DesignTunnel", "Des",
                "Brep av prosjektert tunnel (innvendig geometri, åpen/lukket).",
                GH_ParamAccess.item);

            p.AddCurveParameter("Centerline", "CL",
                "Senterlinje gjennom tunnelen (Curve).",
                GH_ParamAccess.item);

            p.AddPointParameter("A", "A",
                "Startpunkt nær senterlinjen (projiseres på CL).",
                GH_ParamAccess.item);

            p.AddPointParameter("B", "B",
                "Sluttpunkt nær senterlinjen (projiseres på CL).",
                GH_ParamAccess.item);

            p.AddNumberParameter("Step [m]", "dL",
                "Avstand mellom tverrsnitt langs senterlinjen.",
                GH_ParamAccess.item, 1.0);

            p.AddNumberParameter("Tolerance", "Tol",
                "Toleranse (<=0 = dokumenttoleranse).",
                GH_ParamAccess.item, -1.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddNumberParameter("BlastVolume", "V",
                "Tilnærmet volum som må sprenges (sum av alle sprengvolum-solider A–B).",
                GH_ParamAccess.item);

            p.AddBrepParameter("BlastParts", "B",
                "3D-Breps (solider) for hver sammenhengende spreng-sone mellom A og B.",
                GH_ParamAccess.list);

            p.AddBrepParameter("StepBreps", "StepB",
                "Brep per intervall mellom to snitt (step). Kan være null der det ikke er sprengvolum.",
                GH_ParamAccess.list);

            p.AddPointParameter("SectionPoints", "SecPts",
                "Senterlinjepunkt for hvert snitt (XYZ).",
                GH_ParamAccess.list);

            p.AddNumberParameter("SectionBlastArea", "SecA_B",
                "Areal som skal sprenges i hvert tverrsnitt (Design minus Scan).",
                GH_ParamAccess.list);

            p.AddNumberParameter("SectionScanArea", "SecA_S",
                "Areal av scan-profilen i hvert tverrsnitt.",
                GH_ParamAccess.list);

            p.AddNumberParameter("SectionDesignArea", "SecA_D",
                "Areal av design-profilen i hvert tverrsnitt.",
                GH_ParamAccess.list);

            p.AddTextParameter("Info", "i",
                "Diagnose og status.",
                GH_ParamAccess.list);
        }

        // ----------------- INTERN HJELPEKLASSE -----------------

        private class SectionRegion
        {
            public int SectionIndex;
            public double Station;
            public Plane SectionPlane;
            public Curve Curve3D;

            public double MinX, MaxX;
            public double MinY, MaxY;

            public double CenterX, CenterY;
            public double Area;

            public int ClusterId = -1;
        }

        // ----------------- HOVEDLOGIKK -----------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Mesh scanMesh = null;
            Brep designBrep = null;
            Curve centerline = null;
            Point3d ptA = Point3d.Unset;
            Point3d ptB = Point3d.Unset;
            double step = 1.0;
            double tolIn = -1.0;

            if (!da.GetData(0, ref scanMesh) || scanMesh == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler ScanMesh.");
                return;
            }
            if (!da.GetData(1, ref designBrep) || designBrep == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler DesignTunnel.");
                return;
            }
            if (!da.GetData(2, ref centerline) || centerline == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler Centerline.");
                return;
            }
            if (!da.GetData(3, ref ptA))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler punkt A.");
                return;
            }
            if (!da.GetData(4, ref ptB))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler punkt B.");
                return;
            }

            da.GetData(5, ref step);
            da.GetData(6, ref tolIn);

            var info = new List<string>();
            double tol = tolIn <= 0
                ? (RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3)
                : tolIn;
            if (step <= 0) step = 1.0;

            info.Add($"Tol = {tol:0.###}");
            info.Add($"Step = {step:0.###} m");

            // --- Finn A–B langs centerline ---
            if (!centerline.ClosestPoint(ptA, out double tA) ||
                !centerline.ClosestPoint(ptB, out double tB))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Kunne ikke projisere A og/eller B på centerline.");
                return;
            }

            double sA = centerline.GetLength(new Interval(centerline.Domain.Min, tA));
            double sB = centerline.GetLength(new Interval(centerline.Domain.Min, tB));

            double s0 = Math.Min(sA, sB);
            double s1 = Math.Max(sA, sB);

            if (s1 - s0 <= RhinoMath.ZeroTolerance)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A og B er for nære hverandre langs centerline.");
                return;
            }

            double clLen = centerline.GetLength();
            info.Add($"Centerline total lengde: {clLen:0.###} m");
            info.Add($"Analyse mellom stasjon {s0:0.###} m og {s1:0.###} m (A–B).");

            // Stasjoner
            int nSeg = Math.Max(1, (int)Math.Ceiling((s1 - s0) / step));
            var stations = new List<double>();
            for (int i = 0; i <= nSeg; i++)
            {
                double s = s0 + i * step;
                if (s > s1) s = s1;
                stations.Add(s);
            }

            var regionsPerSection = new List<SectionRegion>[stations.Count];

            var sectionPoints = new Point3d[stations.Count];
            var sectionBlastArea = new double[stations.Count];
            var sectionScanArea = new double[stations.Count];
            var sectionDesignArea = new double[stations.Count];

            // NEW: profiler per snitt til viewer
            var scanProfiles2D = new List<Polyline>[stations.Count];
            var designProfiles2D = new List<Polyline>[stations.Count];
            var blastProfiles2D = new List<Polyline>[stations.Count];

            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, Math.Min(4, Environment.ProcessorCount - 1))
            };

            // ---------- PARALLELL SNITT-BEREGNING ----------
            Parallel.For(0, stations.Count, options, i =>
            {
                double s = stations[i];
                var localRegions = new List<SectionRegion>();

                var localScanPolys = new List<Polyline>();
                var localDesPolys = new List<Polyline>();
                var localBlastPolys = new List<Polyline>();

                sectionPoints[i] = Point3d.Unset;
                sectionBlastArea[i] = 0.0;
                sectionScanArea[i] = 0.0;
                sectionDesignArea[i] = 0.0;

                if (!centerline.LengthParameter(s, out double t))
                {
                    regionsPerSection[i] = localRegions;
                    scanProfiles2D[i] = localScanPolys;
                    designProfiles2D[i] = localDesPolys;
                    blastProfiles2D[i] = localBlastPolys;
                    return;
                }

                sectionPoints[i] = centerline.PointAt(t);

                if (!centerline.PerpendicularFrameAt(t, out Plane secPlane))
                {
                    regionsPerSection[i] = localRegions;
                    scanProfiles2D[i] = localScanPolys;
                    designProfiles2D[i] = localDesPolys;
                    blastProfiles2D[i] = localBlastPolys;
                    return;
                }

                // Snitt ScanMesh
                Curve[] scanCrvs = null;
                try { scanCrvs = Mesh.CreateContourCurves(scanMesh, secPlane, tol); }
                catch { scanCrvs = null; }

                // Snitt DesignBrep
                Curve[] desCrvs = null;
                Point3d[] dummyPts;
                try { Intersection.BrepPlane(designBrep, secPlane, tol, out desCrvs, out dummyPts); }
                catch { desCrvs = null; }

                Rhino.Geometry.Transform toXY =
                    Rhino.Geometry.Transform.PlaneToPlane(secPlane, Plane.WorldXY);

                // Areal: scan
                double scanArea = 0.0;
                if (scanCrvs != null)
                {
                    var scanClosed = new List<Curve>();
                    foreach (var c in scanCrvs)
                    {
                        if (c == null || !c.IsClosed) continue;
                        scanClosed.Add(c);

                        var cXY = c.DuplicateCurve();
                        cXY.Transform(toXY);
                        localScanPolys.Add(CurveToPolylineXY(cXY));
                    }

                    if (scanClosed.Count > 0)
                    {
                        var ampS = AreaMassProperties.Compute(scanClosed);
                        if (ampS != null) scanArea = ampS.Area;
                    }
                }

                // Areal: design
                double designArea = 0.0;
                if (desCrvs != null)
                {
                    var desClosed = new List<Curve>();
                    foreach (var c in desCrvs)
                    {
                        if (c == null || !c.IsClosed) continue;
                        desClosed.Add(c);

                        var cXY = c.DuplicateCurve();
                        cXY.Transform(toXY);
                        localDesPolys.Add(CurveToPolylineXY(cXY));
                    }

                    if (desClosed.Count > 0)
                    {
                        var ampD = AreaMassProperties.Compute(desClosed);
                        if (ampD != null) designArea = ampD.Area;
                    }
                }

                // Design MINUS Scan
                var diffCurves = CreatePlanarDifferenceCurves_DesignMinusScan(desCrvs, scanCrvs, tol);

                double blastArea = 0.0;
                if (diffCurves != null && diffCurves.Count > 0)
                {
                    var ampB = AreaMassProperties.Compute(diffCurves);
                    if (ampB != null) blastArea = ampB.Area;

                    foreach (var c in diffCurves)
                    {
                        if (c == null || !c.IsClosed) continue;

                        var cXY = c.DuplicateCurve();
                        cXY.Transform(toXY);
                        localBlastPolys.Add(CurveToPolylineXY(cXY));

                        var bb = cXY.GetBoundingBox(true);
                        if (!bb.IsValid) continue;

                        var ampR = AreaMassProperties.Compute(cXY);
                        if (ampR == null || ampR.Area <= RhinoMath.ZeroTolerance)
                            continue;

                        var r = new SectionRegion
                        {
                            SectionIndex = i,
                            Station = s,
                            SectionPlane = secPlane,
                            Curve3D = c,
                            MinX = bb.Min.X,
                            MaxX = bb.Max.X,
                            MinY = bb.Min.Y,
                            MaxY = bb.Max.Y,
                            CenterX = ampR.Centroid.X,
                            CenterY = ampR.Centroid.Y,
                            Area = ampR.Area
                        };
                        localRegions.Add(r);
                    }
                }

                sectionBlastArea[i] = blastArea;
                sectionScanArea[i] = scanArea;
                sectionDesignArea[i] = designArea;

                regionsPerSection[i] = localRegions;
                scanProfiles2D[i] = localScanPolys;
                designProfiles2D[i] = localDesPolys;
                blastProfiles2D[i] = localBlastPolys;
            });
            // ---------- SLUTT PARALLELL DEL ----------

            var allRegions = new List<SectionRegion>();
            foreach (var lst in regionsPerSection)
                if (lst != null && lst.Count > 0)
                    allRegions.AddRange(lst);

            info.Add($"Antall snitt (A–B): {stations.Count}");
            info.Add($"Snitt med spreng-differanse: {regionsPerSection.Count(r => r != null && r.Count > 0)}");

            if (allRegions.Count == 0)
            {
                info.Add("Fant ingen Design-minus-Scan-områder i tverrsnitt mellom A og B.");

                da.SetData(0, 0.0);
                da.SetDataList(1, new List<Brep>());
                da.SetDataList(2, new List<Brep>());
                da.SetDataList(3, sectionPoints);
                da.SetDataList(4, sectionBlastArea);
                da.SetDataList(5, sectionScanArea);
                da.SetDataList(6, sectionDesignArea);
                da.SetDataList(7, info);
                return;
            }

            // --- Cluster / loft / stepBreps ---

            int clusterCounter = 0;

            if (regionsPerSection[0] != null && regionsPerSection[0].Count > 0)
            {
                foreach (var r in regionsPerSection[0])
                    r.ClusterId = clusterCounter++;
            }

            for (int i = 1; i < regionsPerSection.Length; i++)
            {
                var prev = regionsPerSection[i - 1] ?? new List<SectionRegion>();
                var curr = regionsPerSection[i] ?? new List<SectionRegion>();

                if (curr.Count == 0)
                    continue;

                var prevUsed = new bool[prev.Count];

                foreach (var r in curr)
                {
                    int bestIdx = -1;
                    double bestD2 = double.MaxValue;

                    for (int pi = 0; pi < prev.Count; pi++)
                    {
                        if (prevUsed[pi]) continue;

                        var p = prev[pi];

                        double overlapArea = BBoxOverlapArea(r, p);
                        if (overlapArea <= 0.0) continue;

                        double areaBoxR = (r.MaxX - r.MinX) * (r.MaxY - r.MinY);
                        double areaBoxP = (p.MaxX - p.MinX) * (p.MaxY - p.MinY);
                        double minBox = Math.Max(Math.Min(areaBoxR, areaBoxP), RhinoMath.ZeroTolerance);
                        double frac = overlapArea / minBox;
                        if (frac < 0.2) continue;

                        double dx = r.CenterX - p.CenterX;
                        double dy = r.CenterY - p.CenterY;
                        double d2 = dx * dx + dy * dy;

                        if (d2 < bestD2)
                        {
                            bestD2 = d2;
                            bestIdx = pi;
                        }
                    }

                    if (bestIdx >= 0)
                    {
                        var pBest = prev[bestIdx];

                        double sizeP = Math.Max(pBest.MaxX - pBest.MinX, pBest.MaxY - pBest.MinY);
                        double sizeR = Math.Max(r.MaxX - r.MinX, r.MaxY - r.MinY);
                        double typicalSize = Math.Max(sizeP, sizeR);

                        double maxShift = typicalSize * 0.5;
                        double dist = Math.Sqrt(bestD2);

                        double areaP = Math.Max(pBest.Area, RhinoMath.ZeroTolerance);
                        double areaR = Math.Max(r.Area, RhinoMath.ZeroTolerance);
                        double areaRatio = Math.Min(areaP, areaR) / Math.Max(areaP, areaR);

                        bool okShift = dist <= maxShift;
                        bool okArea = areaRatio >= 0.05;

                        if (okShift && okArea)
                        {
                            r.ClusterId = pBest.ClusterId;
                            prevUsed[bestIdx] = true;
                        }
                        else
                        {
                            r.ClusterId = clusterCounter++;
                        }
                    }
                    else
                    {
                        r.ClusterId = clusterCounter++;
                    }
                }
            }

            var clusters = new Dictionary<int, List<SectionRegion>>();
            foreach (var r in allRegions)
            {
                if (r.ClusterId < 0) continue;
                if (!clusters.TryGetValue(r.ClusterId, out var list))
                {
                    list = new List<SectionRegion>();
                    clusters[r.ClusterId] = list;
                }
                list.Add(r);
            }

            info.Add($"Antall kluster (sammenhengende spreng-soner A–B): {clusters.Count}");

            var blastBreps = new List<Brep>();
            double totalVolume = 0.0;

            int stepCount = Math.Max(0, stations.Count - 1);
            var stepPieces = new List<Brep>[stepCount];
            for (int i = 0; i < stepCount; i++)
                stepPieces[i] = new List<Brep>();

            foreach (var kv in clusters)
            {
                var regs = kv.Value.OrderBy(r => r.SectionIndex).ToList();
                if (regs.Count < 2) continue;

                var pieces = new List<Brep>();

                for (int i = 0; i < regs.Count - 1; i++)
                {
                    var r0 = regs[i];
                    var r1 = regs[i + 1];

                    if (r1.SectionIndex != r0.SectionIndex + 1)
                        continue;

                    var c0 = r0.Curve3D;
                    var c1 = r1.Curve3D;
                    if (c0 == null || c1 == null) continue;

                    if (!Curve.DoDirectionsMatch(c0, c1)) c1.Reverse();

                    var lofts = Brep.CreateFromLoft(
                        new Curve[] { c0, c1 },
                        Point3d.Unset,
                        Point3d.Unset,
                        LoftType.Normal,
                        false);

                    if (lofts == null || lofts.Length == 0) continue;

                    var loftBrep = lofts[0];

                    if (i == 0)
                    {
                        var caps0 = Brep.CreatePlanarBreps(c0, tol);
                        if (caps0 != null && caps0.Length > 0)
                            loftBrep = Brep.JoinBreps(new List<Brep> { loftBrep }.Concat(caps0), tol)?.FirstOrDefault() ?? loftBrep;
                    }
                    if (i == regs.Count - 2)
                    {
                        var caps1 = Brep.CreatePlanarBreps(c1, tol);
                        if (caps1 != null && caps1.Length > 0)
                            loftBrep = Brep.JoinBreps(new List<Brep> { loftBrep }.Concat(caps1), tol)?.FirstOrDefault() ?? loftBrep;
                    }

                    pieces.Add(loftBrep);

                    int segIndex = r0.SectionIndex;
                    if (segIndex >= 0 && segIndex < stepCount)
                        stepPieces[segIndex].Add(loftBrep);
                }

                if (pieces.Count == 0) continue;

                Brep clusterBrep;
                var joined = Brep.JoinBreps(pieces, tol);
                if (joined != null && joined.Length > 0)
                    clusterBrep = joined[0];
                else
                    clusterBrep = pieces[0];

                var cappedCluster = clusterBrep.CapPlanarHoles(tol);
                if (cappedCluster != null && cappedCluster.IsValid)
                    clusterBrep = cappedCluster;

                var vmp = VolumeMassProperties.Compute(clusterBrep);
                if (vmp == null || vmp.Volume <= RhinoMath.ZeroTolerance)
                    continue;

                totalVolume += vmp.Volume;
                blastBreps.Add(clusterBrep);
            }

            var stepBreps = new List<Brep>();
            for (int i = 0; i < stepCount; i++)
            {
                Brep stepB = null;
                var lst = stepPieces[i];
                if (lst != null && lst.Count > 0)
                {
                    var joined = Brep.JoinBreps(lst, tol);
                    stepB = (joined != null && joined.Length > 0) ? joined[0] : lst[0];

                    var capped = stepB.CapPlanarHoles(tol);
                    if (capped != null && capped.IsValid)
                        stepB = capped;
                }
                stepBreps.Add(stepB);
            }

            info.Add($"Antall blast-solider (A–B): {blastBreps.Count}");
            info.Add($"BlastVolume (seksjons-basert, A–B) = {totalVolume:0.###}");

            // lagre til WPF-viewer
            _lastStations = stations;
            _lastSectionPoints = new List<Point3d>(sectionPoints);
            _lastBlastAreas = new List<double>(sectionBlastArea);
            _lastScanAreas = new List<double>(sectionScanArea);
            _lastDesignAreas = new List<double>(sectionDesignArea);
            _lastScanProfiles2D = scanProfiles2D.Select(l => l ?? new List<Polyline>()).ToList();
            _lastDesignProfiles2D = designProfiles2D.Select(l => l ?? new List<Polyline>()).ToList();
            _lastBlastProfiles2D = blastProfiles2D.Select(l => l ?? new List<Polyline>()).ToList();

            // Outputs
            da.SetData(0, totalVolume);
            da.SetDataList(1, blastBreps);
            da.SetDataList(2, stepBreps);
            da.SetDataList(3, sectionPoints);
            da.SetDataList(4, sectionBlastArea);
            da.SetDataList(5, sectionScanArea);
            da.SetDataList(6, sectionDesignArea);
            da.SetDataList(7, info);
        }

        // ---------- HJELP: BBOX-OVERLAPP, POLYLINE, osv. ----------

        private static double BBoxOverlapArea(SectionRegion a, SectionRegion b)
        {
            double x0 = Math.Max(a.MinX, b.MinX);
            double x1 = Math.Min(a.MaxX, b.MaxX);
            double y0 = Math.Max(a.MinY, b.MinY);
            double y1 = Math.Min(a.MaxY, b.MaxY);
            double dx = x1 - x0;
            double dy = y1 - y0;
            if (dx <= 0.0 || dy <= 0.0) return 0.0;
            return dx * dy;
        }

        // NEW: enkel sampling til polyline i XY (for WPF-tegning)
        private static Polyline CurveToPolylineXY(Curve cXY)
        {
            int n = 64;
            if (n < 4) n = 4;

            var pts = new List<Point3d>(n);
            for (int i = 0; i < n; i++)
            {
                double t = cXY.Domain.ParameterAt(i / (double)(n - 1));
                pts.Add(cXY.PointAt(t));
            }
            return new Polyline(pts);
        }

        private static List<Curve> CreatePlanarDifferenceCurves_DesignMinusScan(
            IList<Curve> designCrvs,
            IList<Curve> scanCrvs,
            double tol)
        {
            var result = new List<Curve>();

            var desClosed = new List<Curve>();
            if (designCrvs != null)
            {
                foreach (var c in designCrvs)
                    if (c != null && c.IsClosed)
                        desClosed.Add(c);
            }

            var scanClosed = new List<Curve>();
            if (scanCrvs != null)
            {
                foreach (var c in scanCrvs)
                    if (c != null && c.IsClosed)
                        scanClosed.Add(c);
            }

            if (desClosed.Count == 0)
                return result;

            if (scanClosed.Count == 0)
            {
                result.AddRange(desClosed);
                return result;
            }

            try
            {
                foreach (var d in desClosed)
                {
                    if (d == null) continue;
                    Curve[] diffs = new Curve[] { d };

                    foreach (var s in scanClosed)
                    {
                        if (s == null) continue;
                        var tmp = new List<Curve>();
                        foreach (var c in diffs)
                        {
                            if (c == null) continue;
                            var partial = Curve.CreateBooleanDifference(c, s, tol);
                            if (partial != null && partial.Length > 0)
                                tmp.AddRange(partial);
                        }
                        diffs = tmp.ToArray();
                        if (diffs.Length == 0) break;
                    }

                    foreach (var dd in diffs)
                    {
                        if (dd == null || !dd.IsClosed) continue;
                        var a = AreaMassProperties.Compute(dd);
                        if (a == null || a.Area <= RhinoMath.ZeroTolerance) continue;
                        result.Add(dd);
                    }
                }
            }
            catch
            {
                result.Clear();
                result.AddRange(desClosed);
            }

            return result;
        }

        // ----------------- CONTEXT-MENY / WPF VIEWER -----------------

        protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
        {
            base.AppendAdditionalComponentMenuItems(menu);
            Menu_AppendItem(menu, "Open Section Viewer", OnOpenSectionViewer);
        }

        private void OnOpenSectionViewer(object sender, EventArgs e)
        {
            if (_lastStations == null || _lastStations.Count == 0)
            {
                Rhino.RhinoApp.WriteLine("Ingen seksjonsdata ennå – kjør komponenten først.");
                return;
            }

            var wnd = new SectionViewerWindow(
                _lastStations,
                _lastSectionPoints,
                _lastBlastAreas,
                _lastScanAreas,
                _lastDesignAreas,
                _lastScanProfiles2D,
                _lastDesignProfiles2D,
                _lastBlastProfiles2D);

            wnd.Show();
        }
    }

    // ----------------- WPF-VINDU FOR SEKSJONSVISNING -----------------

    public class SectionViewerWindow : System.Windows.Window
    {
        private readonly List<double> _stations;
        private readonly List<Point3d> _points;
        private readonly List<double> _blastAreas;
        private readonly List<double> _scanAreas;
        private readonly List<double> _designAreas;

        private readonly List<List<Polyline>> _scanProfiles2D;
        private readonly List<List<Polyline>> _designProfiles2D;
        private readonly List<List<Polyline>> _blastProfiles2D;

        private int _index = 0;

        private System.Windows.Controls.TextBox _indexBox;
        private System.Windows.Controls.TextBlock _infoBlock;
        private System.Windows.Controls.TextBlock _blastBlock;

        private System.Windows.Controls.Canvas _canvas; // NEW: tegneflate for tverrsnitt

        public SectionViewerWindow(
            List<double> stations,
            List<Point3d> points,
            List<double> blastAreas,
            List<double> scanAreas,
            List<double> designAreas,
            List<List<Polyline>> scanProfiles2D,
            List<List<Polyline>> designProfiles2D,
            List<List<Polyline>> blastProfiles2D)
        {
            _stations = stations;
            _points = points;
            _blastAreas = blastAreas;
            _scanAreas = scanAreas;
            _designAreas = designAreas;
            _scanProfiles2D = scanProfiles2D;
            _designProfiles2D = designProfiles2D;
            _blastProfiles2D = blastProfiles2D;

            Title = "Tunnel Section Viewer";
            Width = 420;
            Height = 360;

            var root = new System.Windows.Controls.StackPanel
            {
                Margin = new System.Windows.Thickness(10)
            };
            Content = root;

            var topPanel = new System.Windows.Controls.StackPanel
            {
                Orientation = System.Windows.Controls.Orientation.Horizontal,
                Margin = new System.Windows.Thickness(0, 0, 0, 10)
            };
            root.Children.Add(topPanel);

            var prevBtn = new System.Windows.Controls.Button
            {
                Content = "<",
                Width = 40,
                Margin = new System.Windows.Thickness(0, 0, 5, 0)
            };
            prevBtn.Click += (s, e) => { ChangeIndex(_index - 1); };
            topPanel.Children.Add(prevBtn);

            _indexBox = new System.Windows.Controls.TextBox
            {
                Width = 60,
                Margin = new System.Windows.Thickness(0, 0, 5, 0),
                Text = "0"
            };
            _indexBox.KeyDown += (s, e) =>
            {
                if (e.Key == System.Windows.Input.Key.Enter)
                {
                    if (int.TryParse(_indexBox.Text, out int idx))
                        ChangeIndex(idx);
                }
            };
            topPanel.Children.Add(_indexBox);

            var nextBtn = new System.Windows.Controls.Button
            {
                Content = ">",
                Width = 40
            };
            nextBtn.Click += (s, e) => { ChangeIndex(_index + 1); };
            topPanel.Children.Add(nextBtn);

            // NEW: canvas for å tegne tverrsnitt
            _canvas = new System.Windows.Controls.Canvas
            {
                Width = 380,
                Height = 200,
                Background = System.Windows.Media.Brushes.WhiteSmoke,
                Margin = new System.Windows.Thickness(0, 0, 0, 10)
            };
            root.Children.Add(_canvas);

            _infoBlock = new System.Windows.Controls.TextBlock
            {
                TextWrapping = System.Windows.TextWrapping.Wrap,
                Margin = new System.Windows.Thickness(0, 0, 0, 5)
            };
            root.Children.Add(_infoBlock);

            _blastBlock = new System.Windows.Controls.TextBlock
            {
                TextWrapping = System.Windows.TextWrapping.Wrap,
                Foreground = System.Windows.Media.Brushes.Red,
                FontWeight = System.Windows.FontWeights.Bold
            };
            root.Children.Add(_blastBlock);

            ChangeIndex(0);
        }

        private void ChangeIndex(int newIndex)
        {
            if (_stations == null || _stations.Count == 0)
                return;

            if (newIndex < 0) newIndex = 0;
            if (newIndex >= _stations.Count) newIndex = _stations.Count - 1;

            _index = newIndex;
            _indexBox.Text = _index.ToString();
            UpdateInfo();
            DrawSection();   // NEW: tegn valgt snitt
        }

        private void UpdateInfo()
        {
            if (_index < 0 || _index >= _stations.Count)
            {
                _infoBlock.Text = "Index ute av range.";
                _blastBlock.Text = "";
                return;
            }

            double s = _stations[_index];
            Point3d pt = (_points != null && _points.Count > _index)
                ? _points[_index]
                : Point3d.Unset;
            double aB = (_blastAreas != null && _blastAreas.Count > _index) ? _blastAreas[_index] : 0.0;
            double aS = (_scanAreas != null && _scanAreas.Count > _index) ? _scanAreas[_index] : 0.0;
            double aD = (_designAreas != null && _designAreas.Count > _index) ? _designAreas[_index] : 0.0;

            _infoBlock.Text =
                $"Index:        {_index}\n" +
                $"Stasjon:      {s:0.###} m\n" +
                $"CL-punkt:     ({pt.X:0.###}, {pt.Y:0.###}, {pt.Z:0.###})\n" +
                $"A_scan:       {aS:0.###} m²\n" +
                $"A_design:     {aD:0.###} m²";

            _blastBlock.Text =
                $"A_spreng:     {aB:0.###} m²";
        }

        // NEW: tegn tverrsnittet for valgt index
        private void DrawSection()
        {
            if (_canvas == null) return;
            _canvas.Children.Clear();

            if (_index < 0 || _index >= _stations.Count)
                return;

            var scan = (_scanProfiles2D != null && _scanProfiles2D.Count > _index) ? _scanProfiles2D[_index] : null;
            var design = (_designProfiles2D != null && _designProfiles2D.Count > _index) ? _designProfiles2D[_index] : null;
            var blast = (_blastProfiles2D != null && _blastProfiles2D.Count > _index) ? _blastProfiles2D[_index] : null;

            var all = new List<Polyline>();
            if (scan != null) all.AddRange(scan);
            if (design != null) all.AddRange(design);
            if (blast != null) all.AddRange(blast);
            if (all.Count == 0) return;

            double minX = double.PositiveInfinity;
            double maxX = double.NegativeInfinity;
            double minY = double.PositiveInfinity;
            double maxY = double.NegativeInfinity;

            foreach (var pl in all)
            {
                if (pl == null || pl.Count == 0) continue;
                foreach (var p in pl)
                {
                    if (p.X < minX) minX = p.X;
                    if (p.X > maxX) maxX = p.X;
                    if (p.Y < minY) minY = p.Y;
                    if (p.Y > maxY) maxY = p.Y;
                }
            }

            if (!double.IsFinite(minX) || !double.IsFinite(minY)) return;
            if (maxX <= minX || maxY <= minY) return;

            double width = _canvas.Width;
            double height = _canvas.Height;
            double margin = 10.0;

            double sx = (width - 2 * margin) / (maxX - minX);
            double sy = (height - 2 * margin) / (maxY - minY);
            double s = Math.Min(sx, sy);

            Func<Point3d, System.Windows.Point> map = p =>
            {
                double px = (p.X - minX) * s + margin;
                double py = (maxY - p.Y) * s + margin; // inverter Y
                return new System.Windows.Point(px, py);
            };

            void AddPolyList(List<Polyline> polys,
                             System.Windows.Media.Brush stroke,
                             double thickness,
                             System.Windows.Media.Brush fill,
                             bool closed)
            {
                if (polys == null) return;

                foreach (var pl in polys)
                {
                    if (pl == null || pl.Count < 2) continue;

                    var geo = new System.Windows.Media.StreamGeometry();
                    using (var ctx = geo.Open())
                    {
                        var first = map(pl[0]);
                        ctx.BeginFigure(first, fill != null, closed);
                        for (int i = 1; i < pl.Count; i++)
                        {
                            var pt = map(pl[i]);
                            ctx.LineTo(pt, true, false);
                        }
                    }
                    geo.Freeze();

                    var path = new System.Windows.Shapes.Path
                    {
                        Stroke = stroke,
                        StrokeThickness = thickness,
                        Fill = fill,
                        Data = geo
                    };
                    _canvas.Children.Add(path);
                }
            }

            // Tegn i rekkefølge: scan (grå), design (svart), blast (rød med fyll)
            AddPolyList(scan,
                System.Windows.Media.Brushes.Gray,
                1.0,
                null,
                true);

            AddPolyList(design,
                System.Windows.Media.Brushes.Black,
                1.5,
                null,
                true);

            var redFill = new System.Windows.Media.SolidColorBrush(
                System.Windows.Media.Color.FromArgb(80, 255, 0, 0));

            AddPolyList(blast,
                System.Windows.Media.Brushes.Red,
                1.5,
                redFill,
                true);
        }
    }
}
