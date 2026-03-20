using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelAnalys
{
    // =========================================================================
    //  GH_TunnelSectionAnalysis
    //
    //  Snitt-logikk kopiert direkte fra TunnelCompare (Mesh.CreateContourCurves
    //  ? IsClosed-filtrering ? PlaneToPlane-transform ? AreaMassProperties).
    //
    //  Ny logikk: én samling kan ha FLERE mesh (f.eks. bunn + tak).
    //  Kurvene fra alle mesh i samlingen slås sammen med JoinCurves
    //  slik at bunn+tak kan danne én lukket profil.
    //
    //  Faste inputs (0–4):
    //    0  Senterlinje
    //    1  Start [m]   (int)
    //    2  Slutt [m]   (int)
    //    3  Intervall [m]
    //    4  Utmappe     (tekst, valgfri)
    //
    //  Dynamiske par fra indeks 5:
    //    5  MeshSamling1  (Mesh list)
    //    6  Navn1         (string)
    //    7  MeshSamling2  … osv.
    // =========================================================================
    public class GH_TunnelSectionAnalysis : GH_Component, IGH_VariableParameterComponent
    {
        private const int IDX_CL     = 0;
        private const int IDX_START  = 1;
        private const int IDX_END    = 2;
        private const int IDX_STEP   = 3;
        private const int IDX_OUTDIR = 4;
        private const int IDX_DWGON  = 5;
        private const int FIXED_INPUTS = 6;

        private int _pairCount = 1;

        // Sist beregnede data (til viewer)
        private List<double>               _lastStations;
        private List<Point3d>              _lastSectionPoints;
        private List<string>               _lastCollectionNames;
        private List<List<double>>         _lastAreas;
        private List<List<List<Polyline>>> _lastProfiles2D;

        public GH_TunnelSectionAnalysis()
          : base(
              "Tunnel Snittanalyse",
              "SnittAnalys",
              "Beregner tverrsnitt-areal langs senterlinjen for én eller flere mesh-samlinger. " +
              "Flere mesh i én samling (f.eks. bunn og tak) kobles til én lukket profil.",
              "Tunnel",
              "Analysis")
        { }

        public override Guid ComponentGuid =>
            new Guid("A1B2C3D4-E5F6-7890-ABCD-EF1234567890");

        protected override Bitmap Icon => null;
        public override GH_Exposure Exposure => GH_Exposure.primary;

        // =========================================================================
        //  INPUT / OUTPUT
        // =========================================================================

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter("Senterlinje", "CL",
                "Senterlinje gjennom tunnelen.",
                GH_ParamAccess.item);

            p.AddIntegerParameter("Start [m]", "Fra",
                "Stasjon (meter) langs senterlinjen der analysen starter.",
                GH_ParamAccess.item, 0);

            p.AddIntegerParameter("Slutt [m]", "Til",
                "Stasjon (meter) langs senterlinjen der analysen slutter.",
                GH_ParamAccess.item, 100);

            p.AddNumberParameter("Intervall [m]", "dL",
                "Avstand mellom tverrsnitt langs senterlinjen.",
                GH_ParamAccess.item, 5.0);

            p.AddTextParameter("Utmappe", "Dir",
                "Mappe der CSV og DWG lagres. Tom = Desktop\\TunnelAnalys.",
                GH_ParamAccess.item, "");
            p[p.ParamCount - 1].Optional = true;

            p.AddBooleanParameter("Lag DWG", "DWG?",
                "True = eksporter DWG-fil. False = hopp over DWG-eksport.",
                GH_ParamAccess.item, false);

            AddMeshPair(p, 1);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddTextParameter("CsvFil", "CSV",
                "Sti til eksportert CSV-fil.",
                GH_ParamAccess.item);

            p.AddTextParameter("DwgFil", "DWG",
                "Sti til eksportert DWG-fil med tverrprofiler og arealtekst.",
                GH_ParamAccess.item);

            p.AddCurveParameter("DiagramProfiler", "Diag",
                "Alle tverrprofiler som 2D-polylines stablet i rett rutenett.",
                GH_ParamAccess.list);

            p.AddPointParameter("SnittPunkter", "SecPts",
                "Senterlinjepunkt for hvert beregnet tverrsnitt.",
                GH_ParamAccess.list);

            p.AddPointParameter("FraPunkt", "PtFra",
                "Punkt på senterlinjen ved analyser start.",
                GH_ParamAccess.item);

            p.AddPointParameter("TilPunkt", "PtTil",
                "Punkt på senterlinjen ved analyseslutt.",
                GH_ParamAccess.item);

            p.AddTextParameter("Info", "i",
                "Diagnose og status.",
                GH_ParamAccess.list);
        }

        // =========================================================================
        //  IGH_VariableParameterComponent
        // =========================================================================

        public bool CanInsertParameter(GH_ParameterSide side, int index) => false;
        public bool CanRemoveParameter(GH_ParameterSide side, int index) => false;
        public IGH_Param CreateParameter(GH_ParameterSide side, int index) => null;
        public bool DestroyParameter(GH_ParameterSide side, int index) => false;

        public void VariableParameterMaintenance()
        {
            int desired = ComputeDesiredPairCount();
            if (desired == _pairCount) return;

            if (desired > _pairCount)
            {
                while (_pairCount < desired)
                {
                    _pairCount++;
                    int mi = FIXED_INPUTS + (_pairCount - 1) * 2;
                    int ni = mi + 1;

                    Params.RegisterInputParam(new Param_Mesh
                    {
                        Name        = $"MeshSamling{_pairCount}",
                        NickName    = $"M{_pairCount}",
                        Description = $"Mesh-samling #{_pairCount}. Koble inn ett eller flere mesh.",
                        Access      = GH_ParamAccess.list,
                        Optional    = true
                    }, mi);

                    Params.RegisterInputParam(new Param_String
                    {
                        Name        = $"Navn{_pairCount}",
                        NickName    = $"N{_pairCount}",
                        Description = $"Navn på samling #{_pairCount}.",
                        Access      = GH_ParamAccess.item,
                        Optional    = true
                    }, ni);
                }
            }
            else
            {
                while (_pairCount > Math.Max(1, desired))
                {
                    int mi = FIXED_INPUTS + (_pairCount - 1) * 2;
                    int ni = mi + 1;
                    if (ni < Params.Input.Count) Params.UnregisterInputParameter(Params.Input[ni]);
                    if (mi < Params.Input.Count) Params.UnregisterInputParameter(Params.Input[mi]);
                    _pairCount--;
                }
            }
        }

        private int ComputeDesiredPairCount()
        {
            int mi = FIXED_INPUTS + (_pairCount - 1) * 2;
            int ni = mi + 1;
            bool mf = mi < Params.Input.Count && Params.Input[mi].SourceCount > 0;
            bool nf = ni < Params.Input.Count && Params.Input[ni].SourceCount > 0;
            return (mf && nf) ? _pairCount + 1 : _pairCount;
        }

        // =========================================================================
        //  SOLVE
        // =========================================================================

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve  centerline = null;
            int    startM     = 0;
            int    endM       = 100;
            double stepM      = 5.0;
            string outDir     = "";
            bool   lagDwg     = false;

            if (!da.GetData(IDX_CL, ref centerline) || centerline == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler senterlinje.");
                return;
            }
            da.GetData(IDX_START,  ref startM);
            da.GetData(IDX_END,    ref endM);
            da.GetData(IDX_STEP,   ref stepM);
            da.GetData(IDX_OUTDIR, ref outDir);
            da.GetData(IDX_DWGON,  ref lagDwg);

            if (stepM <= 0) stepM = 1.0;
            if (startM > endM) { int tmp = startM; startM = endM; endM = tmp; }

            // Identisk med TunnelCompare: bruk dokumenttoleranse
            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3;

            if (string.IsNullOrWhiteSpace(outDir))
                outDir = Path.Combine(
                    Environment.GetFolderPath(Environment.SpecialFolder.Desktop),
                    "TunnelAnalys");

            // --- Hent dynamiske mesh-samlinger ---
            var collections = new List<(string Name, List<Mesh> Meshes)>();
            for (int pair = 0; pair < _pairCount; pair++)
            {
                int mi = FIXED_INPUTS + pair * 2;
                int ni = mi + 1;
                var meshList = new List<Mesh>();
                string name  = $"Samling{pair + 1}";
                if (mi < Params.Input.Count) da.GetDataList(mi, meshList);
                if (ni < Params.Input.Count) da.GetData(ni, ref name);
                meshList = meshList.Where(m => m != null).ToList();
                if (meshList.Count > 0)
                    collections.Add((name, meshList));
            }

            if (collections.Count == 0)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Ingen mesh-samlinger koblet til.");

            var info = new List<string>();

            // --- Stasjons-grenser (identisk med TunnelCompare) ---
            double clLen = centerline.GetLength();
            double s0    = Math.Max(0, Math.Min(startM, clLen));
            double s1    = Math.Max(0, Math.Min(endM,   clLen));

            info.Add($"Tol = {tol:0.###}");
            info.Add($"Step = {stepM:0.###} m");
            info.Add($"Senterlinje lengde: {clLen:0.###} m");
            info.Add($"Analyse mellom stasjon {s0:0.###} m og {s1:0.###} m");

            if (s1 - s0 <= RhinoMath.ZeroTolerance)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Fra og Til er for nære hverandre.");
                return;
            }

            // Fra/Til-punkter for visuell kontroll
            Point3d ptFra = Point3d.Unset;
            Point3d ptTil = Point3d.Unset;
            if (centerline.LengthParameter(s0, out double tFra)) ptFra = centerline.PointAt(tFra);
            if (centerline.LengthParameter(s1, out double tTil)) ptTil = centerline.PointAt(tTil);

            // --- Stasjoner (identisk med TunnelCompare) ---
            int nSeg = Math.Max(1, (int)Math.Ceiling((s1 - s0) / stepM));
            var stations = new List<double>();
            for (int i = 0; i <= nSeg; i++)
            {
                double s = s0 + i * stepM;
                if (s > s1) s = s1;
                stations.Add(s);
            }

            info.Add($"Antall snitt: {stations.Count}");

            int nSt  = stations.Count;
            int nCol = collections.Count;

            var sectionPoints = new Point3d[nSt];
            var areas         = new double[nCol][];
            var profiles2D    = new List<Polyline>[nCol][];

            for (int c = 0; c < nCol; c++)
            {
                areas[c]      = new double[nSt];
                profiles2D[c] = new List<Polyline>[nSt];
                for (int i = 0; i < nSt; i++)
                    profiles2D[c][i] = new List<Polyline>();
            }

            // Tellevariable for diagnose
            var rawCount    = new int[nCol];
            var closedCount = new int[nCol];
            var hitCount    = new int[nCol];

            var parallelOpts = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, Math.Min(4, Environment.ProcessorCount - 1))
            };

            // ---- PARALLELL SNITT-LØKKE ----
            // Identisk mønster med TunnelCompare sin parallell-løkke.
            Parallel.For(0, nSt, parallelOpts, i =>
            {
                double s = stations[i];

                // -- Identisk med TunnelCompare --
                if (!centerline.LengthParameter(s, out double t)) return;
                sectionPoints[i] = centerline.PointAt(t);
                if (!centerline.PerpendicularFrameAt(t, out Plane secPlane)) return;

                Transform toXY = Transform.PlaneToPlane(secPlane, Plane.WorldXY);

                for (int c = 0; c < nCol; c++)
                {
                    // Kutt ALLE mesh i samlingen – identisk med TunnelCompare sin
                    // Mesh.CreateContourCurves-linie, bare utvidet til en liste.
                    var allCrvs = new List<Curve>();
                    foreach (var mesh in collections[c].Meshes)
                    {
                        if (mesh == null) continue;
                        Curve[] crvs = null;
                        try { crvs = Mesh.CreateContourCurves(mesh, secPlane, tol); }
                        catch { }
                        if (crvs != null)
                            allCrvs.AddRange(crvs.Where(cv => cv != null));
                    }

                    System.Threading.Interlocked.Add(ref rawCount[c], allCrvs.Count);
                    if (allCrvs.Count == 0) continue;

                    // --- NY LOGIKK FOR BUNN+TAK ---
                    // Dersom samlingen har ett mesh ? kurvene er allerede lukkede
                    // (identisk behandling som TunnelCompare).
                    // Dersom samlingen har flere mesh (bunn+tak) ? join kurvene
                    // slik at de to halvdelene kan danne én lukket profil.
                    Curve[] workCrvs;
                    if (collections[c].Meshes.Count > 1)
                    {
                        // Transformer til XY FØR join (samme som TunnelCompare gjør
                        // FØR IsClosed-sjekk)
                        var allXY = new List<Curve>(allCrvs.Count);
                        foreach (var crv in allCrvs)
                        {
                            var cXY = crv.DuplicateCurve();
                            cXY.Transform(toXY);
                            allXY.Add(cXY);
                        }
                        // Join – bunn og tak kan kobles til lukket profil
                        workCrvs = Curve.JoinCurves(allXY, tol) ?? allXY.ToArray();
                    }
                    else
                    {
                        // Ett mesh: transformer til XY (identisk med TunnelCompare)
                        var allXY = new List<Curve>(allCrvs.Count);
                        foreach (var crv in allCrvs)
                        {
                            var cXY = crv.DuplicateCurve();
                            cXY.Transform(toXY);
                            allXY.Add(cXY);
                        }
                        workCrvs = allXY.ToArray();
                    }

                    // -- Identisk med TunnelCompare: filtrer lukkede kurver --
                    var closedCrvs = new List<Curve>();
                    foreach (var crv in workCrvs)
                    {
                        if (crv == null || !crv.IsClosed) continue;
                        closedCrvs.Add(crv);
                    }

                    System.Threading.Interlocked.Add(ref closedCount[c], closedCrvs.Count);
                    if (closedCrvs.Count == 0) continue;

                    // -- Identisk med TunnelCompare: AreaMassProperties --
                    var amp = AreaMassProperties.Compute(closedCrvs);
                    if (amp == null) continue;

                    areas[c][i] = amp.Area;
                    System.Threading.Interlocked.Increment(ref hitCount[c]);

                    // -- Identisk med TunnelCompare: CurveToPolylineXY --
                    var polys = new List<Polyline>(closedCrvs.Count);
                    foreach (var crv in closedCrvs)
                        polys.Add(CurveToPolylineXY(crv));
                    profiles2D[c][i] = polys;
                }
            });
            // ---- SLUTT PARALLELL ----

            for (int c = 0; c < nCol; c++)
                info.Add($"'{collections[c].Name}': råkurver={rawCount[c]}, " +
                         $"lukkede={closedCount[c]}, snitt m/areal={hitCount[c]}/{nSt}");

            // --- Lagre til viewer ---
            _lastStations        = stations;
            _lastSectionPoints   = new List<Point3d>(sectionPoints);
            _lastCollectionNames = collections.Select(x => x.Name).ToList();
            _lastAreas           = areas.Select(a => new List<double>(a)).ToList();
            _lastProfiles2D      = profiles2D
                .Select(col => col.Select(p => p ?? new List<Polyline>()).ToList())
                .ToList();

            // --- Eksporter ---
            string csvPath = ExportCsv(outDir, stations, sectionPoints, collections, areas, info);
            string dwgPath = lagDwg
                ? ExportDwg(outDir, stations, collections, profiles2D, areas, nSt, nCol, info)
                : string.Empty;

            var    diag    = BuildDiagramPolylines(profiles2D, nSt, nCol);

            da.SetData(0, csvPath);
            da.SetData(1, dwgPath);
            da.SetDataList(2, diag);
            da.SetDataList(3, new List<Point3d>(sectionPoints));
            da.SetData(4, ptFra);
            da.SetData(5, ptTil);
            da.SetDataList(6, info);
        }

        // =========================================================================
        //  HJELP – identisk med TunnelCompare
        // =========================================================================

        // Direkte kopi fra TunnelCompare
        private static Polyline CurveToPolylineXY(Curve cXY)
        {
            int n = 64;
            var pts = new List<Point3d>(n);
            for (int i = 0; i < n; i++)
            {
                double t = cXY.Domain.ParameterAt(i / (double)(n - 1));
                pts.Add(cXY.PointAt(t));
            }
            return new Polyline(pts);
        }

        private static void AddMeshPair(GH_InputParamManager p, int number)
        {
            p.AddMeshParameter(
                $"MeshSamling{number}", $"M{number}",
                $"Mesh-samling #{number}. Koble inn ett eller flere mesh (f.eks. bunn og tak). " +
                "Kurver fra alle mesh joinæs til én lukket profil.",
                GH_ParamAccess.list);
            p[p.ParamCount - 1].Optional = true;

            p.AddTextParameter(
                $"Navn{number}", $"N{number}",
                $"Navn på mesh-samling #{number}.",
                GH_ParamAccess.item, $"Samling{number}");
            p[p.ParamCount - 1].Optional = true;
        }

        // =========================================================================
        //  CSV-EKSPORT
        // =========================================================================

        private static string ExportCsv(
            string                                  outDir,
            List<double>                            stations,
            Point3d[]                               points,
            List<(string Name, List<Mesh> Meshes)> collections,
            double[][]                              areas,
            List<string>                            info)
        {
            try
            {
                Directory.CreateDirectory(outDir);
                string path = Path.Combine(outDir,
                    $"TunnelSnittAnalys_{DateTime.Now:yyyyMMdd_HHmmss}.csv");

                using var sw = new StreamWriter(path, false, System.Text.Encoding.UTF8);

                var header = new List<string> { "Stasjon [m]", "CL-X", "CL-Y", "CL-Z" };
                foreach (var col in collections)
                    header.Add($"Areal_{col.Name} [m2]");
                sw.WriteLine(string.Join(";", header));

                for (int i = 0; i < stations.Count; i++)
                {
                    var row = new List<string>
                    {
                        stations[i].ToString("0.###"),
                        points[i].X.ToString("0.###"),
                        points[i].Y.ToString("0.###"),
                        points[i].Z.ToString("0.###")
                    };
                    for (int c = 0; c < collections.Count; c++)
                        row.Add(areas[c][i].ToString("0.####"));
                    sw.WriteLine(string.Join(";", row));
                }

                info.Add($"CSV: {path}");
                return path;
            }
            catch (Exception ex)
            {
                info.Add($"CSV feil: {ex.Message}");
                return string.Empty;
            }
        }

        // =========================================================================
        //  DWG-EKSPORT
        // =========================================================================

        private static string ExportDwg(
            string                                  outDir,
            List<double>                            stations,
            List<(string Name, List<Mesh> Meshes)> collections,
            List<Polyline>[][]                      profiles2D,
            double[][]                              areas,
            int                                     nSt,
            int                                     nCol,
            List<string>                            info)
        {
            if (nSt == 0 || nCol == 0) return string.Empty;

            try
            {
                double maxSpan = 0.0;
                for (int c = 0; c < nCol; c++)
                    for (int i = 0; i < nSt; i++)
                        foreach (var pl in profiles2D[c][i])
                        {
                            if (pl == null) continue;
                            var bb  = pl.BoundingBox;
                            double sp = Math.Max(bb.Max.X - bb.Min.X, bb.Max.Y - bb.Min.Y);
                            if (sp > maxSpan) maxSpan = sp;
                        }

                if (maxSpan <= RhinoMath.ZeroTolerance) maxSpan = 1.0;

                double cellSize = maxSpan * 1.15;
                double gap      = cellSize * 0.40;
                double colStep  = cellSize + gap;
                double rowStep  = cellSize + gap;
                double textH    = Math.Max(cellSize * 0.07, 0.05);

                var tempDoc = RhinoDoc.CreateHeadless(null);
                if (tempDoc == null) { info.Add("DWG: kunne ikke opprette dokument."); return string.Empty; }

                try
                {
                    for (int c = 0; c < nCol; c++)
                    {
                        double baseY = -c * rowStep;

                        // Samlingsnavnlabel
                        string rowLabel = c < collections.Count ? collections[c].Name : $"Samling{c + 1}";
                        tempDoc.Objects.AddText(new Rhino.Geometry.TextEntity
                        {
                            PlainText     = rowLabel,
                            Plane         = new Plane(new Point3d(-gap * 0.8, baseY + cellSize * 0.5, 0), Vector3d.ZAxis),
                            TextHeight    = textH,
                            Justification = TextJustification.MiddleRight
                        });

                        for (int i = 0; i < nSt; i++)
                        {
                            double baseX = i * colStep;

                            // Stasjonslabel
                            if (c == 0)
                                tempDoc.Objects.AddText(new Rhino.Geometry.TextEntity
                                {
                                    PlainText     = $"{stations[i]:0.#} m",
                                    Plane         = new Plane(new Point3d(baseX + cellSize * 0.5, gap * 0.75, 0), Vector3d.ZAxis),
                                    TextHeight    = textH,
                                    Justification = TextJustification.BottomCenter
                                });

                            // Celle-ramme
                            var frameCrv = new Polyline(new[]
                            {
                                new Point3d(baseX,            baseY,            0),
                                new Point3d(baseX + cellSize, baseY,            0),
                                new Point3d(baseX + cellSize, baseY + cellSize, 0),
                                new Point3d(baseX,            baseY + cellSize, 0),
                                new Point3d(baseX,            baseY,            0)
                            }).ToNurbsCurve();
                            if (frameCrv != null)
                            {
                                var fa = new Rhino.DocObjects.ObjectAttributes();
                                fa.ColorSource = Rhino.DocObjects.ObjectColorSource.ColorFromObject;
                                fa.ObjectColor = Color.FromArgb(200, 200, 200);
                                tempDoc.Objects.Add(frameCrv, fa);
                            }

                            var polys = profiles2D[c][i];
                            if (polys == null || polys.Count == 0) continue;

                            // Sentrering i cellen
                            double pMinX =  double.PositiveInfinity, pMaxX = double.NegativeInfinity;
                            double pMinY =  double.PositiveInfinity, pMaxY = double.NegativeInfinity;
                            foreach (var pl in polys)
                                if (pl != null)
                                    foreach (var pt in pl)
                                    {
                                        if (pt.X < pMinX) pMinX = pt.X;
                                        if (pt.X > pMaxX) pMaxX = pt.X;
                                        if (pt.Y < pMinY) pMinY = pt.Y;
                                        if (pt.Y > pMaxY) pMaxY = pt.Y;
                                    }
                            if (!double.IsFinite(pMinX)) continue;

                            double offsetX = (baseX + cellSize * 0.5) - (pMinX + pMaxX) * 0.5;
                            double offsetY = (baseY + cellSize * 0.55) - (pMinY + pMaxY) * 0.5;

                            var profAttr = new Rhino.DocObjects.ObjectAttributes();
                            profAttr.ColorSource = Rhino.DocObjects.ObjectColorSource.ColorFromObject;
                            profAttr.ObjectColor = CollectionColor(c);

                            foreach (var pl in polys)
                            {
                                if (pl == null || pl.Count < 2) continue;
                                var moved = new Polyline(pl.Select(pt =>
                                    new Point3d(pt.X + offsetX, pt.Y + offsetY, 0)));
                                var crv = moved.ToNurbsCurve();
                                if (crv != null) tempDoc.Objects.Add(crv, profAttr);
                            }

                            // Arealtekst
                            double areal   = areas[c][i];
                            string areaStr = areal > RhinoMath.ZeroTolerance
                                ? $"A={areal:0.##} m²" : "A=0";
                            tempDoc.Objects.AddText(new Rhino.Geometry.TextEntity
                            {
                                PlainText     = areaStr,
                                Plane         = new Plane(new Point3d(baseX + cellSize * 0.5, baseY + textH * 0.6, 0), Vector3d.ZAxis),
                                TextHeight    = textH * 0.85,
                                Justification = TextJustification.BottomCenter
                            });
                        }
                    }

                    Directory.CreateDirectory(outDir);
                    string path = Path.Combine(outDir,
                        $"TunnelSnittAnalys_{DateTime.Now:yyyyMMdd_HHmmss}.dwg");

                    bool ok = tempDoc.Export(path);
                    if (ok) { info.Add($"DWG: {path}"); return path; }
                    info.Add("DWG: eksport feilet.");
                    return string.Empty;
                }
                finally { tempDoc.Dispose(); }
            }
            catch (Exception ex)
            {
                info.Add($"DWG feil: {ex.Message}");
                return string.Empty;
            }
        }

        private static Color CollectionColor(int index) =>
            (index % 6) switch
            {
                0 => Color.FromArgb(70,  130, 180),
                1 => Color.FromArgb(220, 60,  0),
                2 => Color.FromArgb(46,  139, 87),
                3 => Color.FromArgb(200, 150, 20),
                4 => Color.FromArgb(140, 100, 210),
                _ => Color.FromArgb(0,   128, 128),
            };

        // =========================================================================
        //  DIAGRAM – GH-geometri
        // =========================================================================

        private static List<Curve> BuildDiagramPolylines(
            List<Polyline>[][] profiles2D, int nSt, int nCol)
        {
            var result = new List<Curve>();
            if (nSt == 0 || nCol == 0) return result;

            double maxSpan = 0.0;
            for (int c = 0; c < nCol; c++)
                for (int i = 0; i < nSt; i++)
                    foreach (var pl in profiles2D[c][i])
                    {
                        if (pl == null) continue;
                        var bb = pl.BoundingBox;
                        double sp = Math.Max(bb.Max.X - bb.Min.X, bb.Max.Y - bb.Min.Y);
                        if (sp > maxSpan) maxSpan = sp;
                    }

            if (maxSpan <= RhinoMath.ZeroTolerance) maxSpan = 1.0;
            double cell = maxSpan * 1.3;

            for (int c = 0; c < nCol; c++)
            {
                double baseY = -c * cell;
                for (int i = 0; i < nSt; i++)
                {
                    double baseX = i * cell;
                    foreach (var pl in profiles2D[c][i])
                    {
                        if (pl == null || pl.Count < 2) continue;
                        var nc = new Polyline(
                            pl.Select(p => new Point3d(baseX + p.X, baseY + p.Y, 0)))
                            .ToNurbsCurve();
                        if (nc != null) result.Add(nc);
                    }
                }
            }

            return result;
        }

        // =========================================================================
        //  CONTEXT-MENY
        // =========================================================================

        protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
        {
            base.AppendAdditionalComponentMenuItems(menu);
            Menu_AppendItem(menu, "Åpne Tverrprofil-viewer", OnOpenViewer);
        }

        private void OnOpenViewer(object sender, EventArgs e)
        {
            if (_lastStations == null || _lastStations.Count == 0)
            {
                RhinoApp.WriteLine("Ingen seksjonsdata ennå – kjør komponenten først.");
                return;
            }

            new SectionAnalysisViewerWindow(
                _lastStations,
                _lastSectionPoints,
                _lastCollectionNames,
                _lastAreas,
                _lastProfiles2D).Show();
        }
    }

    // =========================================================================
    //  WPF VIEWER – identisk struktur med TunnelCompare sin SectionViewerWindow
    // =========================================================================

    public class SectionAnalysisViewerWindow : System.Windows.Window
    {
        private readonly List<double>               _stations;
        private readonly List<Point3d>              _points;
        private readonly List<string>               _names;
        private readonly List<List<double>>         _areas;
        private readonly List<List<List<Polyline>>> _profiles2D;

        private int _index = 0;

        private System.Windows.Controls.TextBox   _indexBox;
        private System.Windows.Controls.TextBlock _stationBlock;
        private System.Windows.Controls.TextBlock _areasBlock;
        private System.Windows.Controls.Canvas    _canvas;

        private static readonly System.Windows.Media.Color[] Palette =
        {
            System.Windows.Media.Colors.SteelBlue,
            System.Windows.Media.Colors.OrangeRed,
            System.Windows.Media.Colors.SeaGreen,
            System.Windows.Media.Colors.Goldenrod,
            System.Windows.Media.Colors.MediumPurple,
            System.Windows.Media.Colors.Teal,
        };

        public SectionAnalysisViewerWindow(
            List<double>               stations,
            List<Point3d>              points,
            List<string>               names,
            List<List<double>>         areas,
            List<List<List<Polyline>>> profiles2D)
        {
            _stations   = stations;
            _points     = points;
            _names      = names;
            _areas      = areas;
            _profiles2D = profiles2D;

            Title      = "Tverrprofil-viewer";
            Width      = 580;
            Height     = 580;
            ResizeMode = System.Windows.ResizeMode.CanResize;

            BuildUI();
            ChangeIndex(0);
        }

        private void BuildUI()
        {
            var root = new System.Windows.Controls.Grid();
            root.RowDefinitions.Add(new System.Windows.Controls.RowDefinition { Height = System.Windows.GridLength.Auto });
            root.RowDefinitions.Add(new System.Windows.Controls.RowDefinition { Height = new System.Windows.GridLength(1, System.Windows.GridUnitType.Star) });
            root.RowDefinitions.Add(new System.Windows.Controls.RowDefinition { Height = System.Windows.GridLength.Auto });
            root.RowDefinitions.Add(new System.Windows.Controls.RowDefinition { Height = System.Windows.GridLength.Auto });
            Content = root;

            // Navigasjon
            var nav = new System.Windows.Controls.StackPanel
            {
                Orientation = System.Windows.Controls.Orientation.Horizontal,
                Margin      = new System.Windows.Thickness(8, 8, 8, 4)
            };
            System.Windows.Controls.Grid.SetRow(nav, 0);
            root.Children.Add(nav);

            var btnPrev = new System.Windows.Controls.Button { Content = "?", Width = 36, Height = 28, Margin = new System.Windows.Thickness(0, 0, 4, 0) };
            btnPrev.Click += (s, e) => ChangeIndex(_index - 1);
            nav.Children.Add(btnPrev);

            _indexBox = new System.Windows.Controls.TextBox
            {
                Width = 56, Height = 28,
                VerticalContentAlignment = System.Windows.VerticalAlignment.Center,
                Margin = new System.Windows.Thickness(0, 0, 4, 0),
                Text = "0"
            };
            _indexBox.KeyDown += (s, e) =>
            {
                if (e.Key == System.Windows.Input.Key.Enter &&
                    int.TryParse(_indexBox.Text, out int idx))
                    ChangeIndex(idx);
            };
            nav.Children.Add(_indexBox);

            var btnNext = new System.Windows.Controls.Button { Content = "?", Width = 36, Height = 28, Margin = new System.Windows.Thickness(0, 0, 12, 0) };
            btnNext.Click += (s, e) => ChangeIndex(_index + 1);
            nav.Children.Add(btnNext);

            _stationBlock = new System.Windows.Controls.TextBlock
            {
                VerticalAlignment = System.Windows.VerticalAlignment.Center,
                FontWeight        = System.Windows.FontWeights.SemiBold
            };
            nav.Children.Add(_stationBlock);

            // Canvas
            var border = new System.Windows.Controls.Border
            {
                BorderBrush     = System.Windows.Media.Brushes.LightGray,
                BorderThickness = new System.Windows.Thickness(1),
                Margin          = new System.Windows.Thickness(8, 0, 8, 4),
                Background      = System.Windows.Media.Brushes.WhiteSmoke
            };
            System.Windows.Controls.Grid.SetRow(border, 1);
            root.Children.Add(border);

            _canvas = new System.Windows.Controls.Canvas
            {
                HorizontalAlignment = System.Windows.HorizontalAlignment.Stretch,
                VerticalAlignment   = System.Windows.VerticalAlignment.Stretch
            };
            border.Child = _canvas;
            _canvas.SizeChanged += (s, e) => DrawSection();

            // Areal-tekst
            _areasBlock = new System.Windows.Controls.TextBlock
            {
                Margin       = new System.Windows.Thickness(8, 2, 8, 4),
                TextWrapping = System.Windows.TextWrapping.Wrap,
                FontFamily   = new System.Windows.Media.FontFamily("Consolas")
            };
            System.Windows.Controls.Grid.SetRow(_areasBlock, 2);
            root.Children.Add(_areasBlock);

            // Fargeforklaring
            var legend = new System.Windows.Controls.WrapPanel { Margin = new System.Windows.Thickness(8, 0, 8, 8) };
            System.Windows.Controls.Grid.SetRow(legend, 3);
            root.Children.Add(legend);

            if (_names != null)
                for (int c = 0; c < _names.Count; c++)
                {
                    var col = Palette[c % Palette.Length];
                    legend.Children.Add(new System.Windows.Shapes.Rectangle
                    {
                        Width  = 12, Height = 12,
                        Fill   = new System.Windows.Media.SolidColorBrush(col),
                        Margin = new System.Windows.Thickness(0, 2, 4, 0),
                        VerticalAlignment = System.Windows.VerticalAlignment.Center
                    });
                    legend.Children.Add(new System.Windows.Controls.TextBlock
                    {
                        Text   = _names[c],
                        Margin = new System.Windows.Thickness(0, 2, 14, 0),
                        VerticalAlignment = System.Windows.VerticalAlignment.Center
                    });
                }
        }

        private void ChangeIndex(int newIndex)
        {
            if (_stations == null || _stations.Count == 0) return;
            newIndex       = Math.Max(0, Math.Min(newIndex, _stations.Count - 1));
            _index         = newIndex;
            _indexBox.Text = _index.ToString();
            UpdateInfo();
            DrawSection();
        }

        private void UpdateInfo()
        {
            int cnt = _stations?.Count ?? 0;
            if (_index < 0 || _index >= cnt) { _stationBlock.Text = ""; _areasBlock.Text = ""; return; }

            double sv = _stations[_index];
            var    pt = (_points != null && _points.Count > _index) ? _points[_index] : Point3d.Unset;

            _stationBlock.Text = $"Snitt {_index} / {cnt - 1}  |  {sv:0.###} m  |  ({pt.X:0.##}, {pt.Y:0.##}, {pt.Z:0.##})";

            var sb = new System.Text.StringBuilder();
            if (_names != null)
                for (int c = 0; c < _names.Count; c++)
                {
                    double a = (_areas != null && c < _areas.Count && _index < _areas[c].Count)
                        ? _areas[c][_index] : 0.0;
                    sb.AppendLine($"{_names[c],-20}  {a:0.####} m²");
                }
            _areasBlock.Text = sb.ToString();
        }

        // DrawSection: identisk mønster med TunnelCompare sin DrawSection
        private void DrawSection()
        {
            if (_canvas == null) return;
            _canvas.Children.Clear();

            int cnt = _stations?.Count ?? 0;
            if (_index < 0 || _index >= cnt) return;

            double w = _canvas.ActualWidth;
            double h = _canvas.ActualHeight;
            if (w < 10 || h < 10) return;

            int nCol = _profiles2D?.Count ?? 0;
            var all  = new List<Polyline>();
            for (int c = 0; c < nCol; c++)
            {
                if (_profiles2D[c] == null || _index >= _profiles2D[c].Count) continue;
                var lst = _profiles2D[c][_index];
                if (lst != null) all.AddRange(lst);
            }

            if (all.Count == 0)
            {
                var tb = new System.Windows.Controls.TextBlock
                    { Text = "Ingen profiler for dette snittet.", Foreground = System.Windows.Media.Brushes.Gray };
                System.Windows.Controls.Canvas.SetLeft(tb, w / 2 - 100);
                System.Windows.Controls.Canvas.SetTop(tb,  h / 2 - 10);
                _canvas.Children.Add(tb);
                return;
            }

            double minX =  double.PositiveInfinity, maxX = double.NegativeInfinity;
            double minY =  double.PositiveInfinity, maxY = double.NegativeInfinity;
            foreach (var pl in all)
            {
                if (pl == null) continue;
                foreach (var p in pl)
                {
                    if (p.X < minX) minX = p.X;
                    if (p.X > maxX) maxX = p.X;
                    if (p.Y < minY) minY = p.Y;
                    if (p.Y > maxY) maxY = p.Y;
                }
            }

            if (!double.IsFinite(minX) || maxX <= minX || maxY <= minY) return;

            const double mg = 16.0;
            double sc = Math.Min((w - 2 * mg) / (maxX - minX), (h - 2 * mg) / (maxY - minY));

            // Identisk med TunnelCompare sin map-funksjon
            System.Windows.Point Map(Point3d p) =>
                new System.Windows.Point((p.X - minX) * sc + mg, (maxY - p.Y) * sc + mg);

            // Tegn profiler
            for (int c = 0; c < nCol; c++)
            {
                if (_profiles2D[c] == null || _index >= _profiles2D[c].Count) continue;
                var lst = _profiles2D[c][_index];
                if (lst == null) continue;

                var col    = Palette[c % Palette.Length];
                var stroke = new System.Windows.Media.SolidColorBrush(col);
                var fill   = new System.Windows.Media.SolidColorBrush(
                    System.Windows.Media.Color.FromArgb(40, col.R, col.G, col.B));

                foreach (var pl in lst)
                {
                    if (pl == null || pl.Count < 2) continue;
                    var geo = new System.Windows.Media.StreamGeometry();
                    using (var ctx = geo.Open())
                    {
                        ctx.BeginFigure(Map(pl[0]), true, true);
                        for (int j = 1; j < pl.Count; j++)
                            ctx.LineTo(Map(pl[j]), true, false);
                    }
                    geo.Freeze();
                    _canvas.Children.Add(new System.Windows.Shapes.Path
                    {
                        Stroke = stroke, StrokeThickness = 1.5, Fill = fill, Data = geo
                    });
                }
            }
        }
    }
}
