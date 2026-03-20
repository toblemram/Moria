using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    /// <summary>
    /// SOSI.ImportX (3 moduser i første omgang)
    ///
    /// Mode:
    /// 0 = Generell
    ///     Outputs:
    ///       - Punkt  : alle .PUNKT (raw)
    ///       - Kurver : alle .KURVE (raw polyline)
    ///       - Flater : alle .FLATE (planar Brep hvis mulig)
    ///       - + egne output-lister per OBJTYPE (inneholder alle objekter: punkt/kurver/flater)
    ///         Manglende OBJTYPE -> "Udefinert"
    ///       - Info
    ///
    /// 1 = Rør
    ///     Inputs: Diameter
    ///     Outputs:
    ///       - Rør (Brep) laget fra alle KURVE
    ///       - ObjType (tekst) 1-til-1 med rørene, fallback "Udefinert"
    ///       - Info
    ///
    /// 2 = Boks
    ///     Inputs: Høyde, Bredde, Lengde
    ///     Outputs:
    ///       - Bokser (Brep) laget fra alle PUNKT
    ///       - ObjType (tekst) 1-til-1 med boksene, fallback "Udefinert"
    ///       - Info
    ///
    /// NB: Generell har dynamiske outputs for OBJTYPE. Det betyr at outputs kan endre seg når filen endrer seg.
    /// </summary>
    public class GH_SosiImportX : GH_Component, IGH_VariableParameterComponent
    {
        public GH_SosiImportX() : base(
            "SOSI.ImportX",
            "SOSI.X",
            "SOSI-import med moduser: Generell, Rør, Boks.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("C21B45A5-0D85-47D5-9A21-6BC2C9E2A2F1");
        protected override System.Drawing.Bitmap Icon => null;

        private enum ImportMode
        {
            Generell = 0,
            Ror = 1,
            Boks = 2
        }

        private ImportMode _lastMode = (ImportMode)(-1);
        private string _lastGenerellTypeSig = "";

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);

            p.AddIntegerParameter("Mode", "mode", "0=Generell, 1=Rør, 2=Boks", GH_ParamAccess.item, 0);
            var modeParam = p[1] as Param_Integer;
            if (modeParam != null)
            {
                modeParam.AddNamedValue("Generell", (int)ImportMode.Generell);
                modeParam.AddNamedValue("Rør", (int)ImportMode.Ror);
                modeParam.AddNamedValue("Boks", (int)ImportMode.Boks);
            }

            // Start med Generell inputs (ingen ekstra)
            EnsureInputsForMode(ImportMode.Generell);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            // Start med Generell outputs (minimal sett - blir justert i Solve når vi har types)
            EnsureOutputsForMode(ImportMode.Generell, new List<string> { "Udefinert" });
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string path = null;
            int modeInt = 0;

            if (!da.GetData(0, ref path) || string.IsNullOrWhiteSpace(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler sti til SOSI-fil.");
                return;
            }
            da.GetData(1, ref modeInt);

            var mode = (ImportMode)Math.Max(0, Math.Min(2, modeInt));

            // Mode bytte -> rebuild inputs/outputs og re-solve
            if (mode != _lastMode)
            {
                _lastMode = mode;
                EnsureInputsForMode(mode);
                EnsureOutputsForMode(mode, new List<string> { "Udefinert" });
                Params.OnParametersChanged();
                ExpireSolution(true);
                return;
            }

            if (!File.Exists(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Finner ikke fil: {path}");
                return;
            }

            string[] lines;
            try { lines = File.ReadAllLines(path, Encoding.UTF8); }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Kunne ikke lese fil: {ex.Message}");
                return;
            }

            var data = SosiParser.Parse(lines);

            // Generell: må kanskje oppdatere dynamiske outputs basert på faktisk OBJTYPE i filen
            if (mode == ImportMode.Generell)
            {
                var typeKeys = CollectAllObjTypes(data); // inkluderer Udefinert hvis mangler/noe
                var sig = string.Join("|", typeKeys.Select(t => t.ToLowerInvariant()));

                if (!string.Equals(sig, _lastGenerellTypeSig, StringComparison.Ordinal))
                {
                    _lastGenerellTypeSig = sig;
                    EnsureOutputsForMode(mode, typeKeys);
                    Params.OnParametersChanged();
                    ExpireSolution(true);
                    return;
                }

                SolveGenerell(da, data, typeKeys);
                return;
            }

            if (mode == ImportMode.Ror)
            {
                SolveRor(da, data);
                return;
            }

            if (mode == ImportMode.Boks)
            {
                SolveBoks(da, data);
                return;
            }
        }

        // ------------------------ SOLVE: GENERELL ------------------------

        private void SolveGenerell(IGH_DataAccess da, SosiData d, List<string> typeKeys)
        {
            // Outputs:
            // 0 Punkt (Point list)
            // 1 Kurver (Curve list)
            // 2 Flater (Brep list)
            // 3..(3+typeKeys-1) OBJTYPE buckets (Geometry list)
            // last: Info (Text list)

            var pts = d.Points.Where(p => p.Pos.IsValid).Select(p => p.Pos).ToList();

            var curves = d.Curves
                .Where(c => c.Points.Count >= 2)
                .Select(c => (Curve)new PolylineCurve(c.Points))
                .ToList();

            int flateFail = 0;
            var flater = new List<Brep>();
            var flateBrepsByIndex = new Dictionary<int, Brep>(); // for bucket-add

            for (int i = 0; i < d.Flater.Count; i++)
            {
                var f = d.Flater[i];
                if (f.Points.Count < 3) continue;

                var brep = TryCreatePlanarBrep(f.Points);
                if (brep != null)
                {
                    flateBrepsByIndex[i] = brep;
                    flater.Add(brep);
                }
                else
                {
                    flateFail++;
                }
            }

            // Bygg buckets per OBJTYPE (alle objekter: punkt/kurve/flate)
            // Bucket type: List<object> (Point3d/Curve/Brep) - Param_Geometry klarer mixed.
            var buckets = typeKeys.ToDictionary(k => k, _ => new List<object>(), StringComparer.OrdinalIgnoreCase);

            // Punkt til buckets
            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                string t = NormalizeType(p.ObjType);
                if (!buckets.TryGetValue(t, out var list))
                    list = buckets[t] = new List<object>();
                list.Add(p.Pos);
            }

            // Kurver til buckets (bruk polylinekurve)
            foreach (var c in d.Curves.Where(c => c.Points.Count >= 2))
            {
                string t = NormalizeType(c.ObjType);
                if (!buckets.TryGetValue(t, out var list))
                    list = buckets[t] = new List<object>();
                list.Add(new PolylineCurve(c.Points));
            }

            // Flater til buckets (kun de vi klarte å lage brep for)
            for (int i = 0; i < d.Flater.Count; i++)
            {
                if (!flateBrepsByIndex.TryGetValue(i, out var brep)) continue;

                string t = NormalizeType(d.Flater[i].ObjType);
                if (!buckets.TryGetValue(t, out var list))
                    list = buckets[t] = new List<object>();
                list.Add(brep);
            }

            var info = new List<string>
            {
                "Mode=Generell",
                $"PUNKT={pts.Count}",
                $"KURVE={curves.Count}",
                $".FLATE (brep)={flater.Count} (feilet={flateFail})",
                $"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})",
                $"OBJTYPE outputs={typeKeys.Count} (inkl. 'Udefinert')"
            };

            da.SetDataList(0, pts);
            da.SetDataList(1, curves);
            da.SetDataList(2, flater);

            int start = 3;
            for (int i = 0; i < typeKeys.Count; i++)
            {
                var key = typeKeys[i];
                var list = buckets.TryGetValue(key, out var objs) ? objs : new List<object>();
                da.SetDataList(start + i, list);
            }

            da.SetDataList(start + typeKeys.Count, info);
        }

        // ------------------------ SOLVE: RØR ------------------------

        private void SolveRor(IGH_DataAccess da, SosiData d)
        {
            // Inputs: [0]=path, [1]=mode, [2]=diameter
            double diameter = 1.0;
            da.GetData(2, ref diameter);
            double radius = Math.Max(1e-6, diameter * 0.5);

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;

            var pipes = new List<Brep>();
            var types = new List<string>();

            int fail = 0;

            foreach (var c in d.Curves.Where(c => c.Points.Count >= 2))
            {
                var path = (Curve)new PolylineCurve(c.Points);

                var pp = Brep.CreatePipe(path, radius, false, PipeCapMode.Round, true, tol, RhinoMath.ToRadians(1.0));
                if (pp != null && pp.Length > 0 && pp[0] != null)
                {
                    pipes.Add(pp[0]);
                    types.Add(NormalizeType(c.ObjType)); // 1-til-1
                }
                else
                {
                    fail++;
                }
            }

            var info = new List<string>
            {
                "Mode=Rør",
                $"KURVE (input)={d.Curves.Count(c => c.Points.Count>=2)}",
                $"Rør (brep)={pipes.Count} (feilet={fail})",
                $"Diameter={diameter}",
                $"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})"
            };

            da.SetDataList(0, pipes);
            da.SetDataList(1, types);
            da.SetDataList(2, info);
        }

        // ------------------------ SOLVE: BOKS ------------------------

        private void SolveBoks(IGH_DataAccess da, SosiData d)
        {
            // Inputs: [0]=path, [1]=mode, [2]=H, [3]=B, [4]=L
            double height = 2.0, width = 1.0, length = 1.0;
            da.GetData(2, ref height);
            da.GetData(3, ref width);
            da.GetData(4, ref length);

            var boxes = new List<Brep>();
            var types = new List<string>();

            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                var pt = p.Pos;
                var plane = new Plane(pt, Vector3d.ZAxis);

                var ix = new Interval(-width * 0.5, width * 0.5);
                var iy = new Interval(-length * 0.5, length * 0.5);
                var iz = new Interval(-height * 0.5, height * 0.5);

                var box = new Box(plane, ix, iy, iz).ToBrep();
                if (box != null)
                {
                    boxes.Add(box);
                    types.Add(NormalizeType(p.ObjType)); // 1-til-1
                }
            }

            var info = new List<string>
            {
                "Mode=Boks",
                $"PUNKT={d.Points.Count(p => p.Pos.IsValid)}",
                $"Bokser={boxes.Count}",
                $"Dim(H,B,L)=({height},{width},{length})",
                $"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})"
            };

            da.SetDataList(0, boxes);
            da.SetDataList(1, types);
            da.SetDataList(2, info);
        }

        // ------------------------ INPUT/OUTPUT SETUP ------------------------

        private void EnsureInputsForMode(ImportMode mode)
        {
            // Hold alltid [0]=path og [1]=mode
            while (Params.Input.Count > 2)
                Params.UnregisterInputParameter(Params.Input[2], true);

            switch (mode)
            {
                case ImportMode.Generell:
                    // ingen ekstra inputs
                    break;

                case ImportMode.Ror:
                    AddNumberInput("Diameter", "D", "Rørdiameter for alle KURVE.", 1.0);
                    break;

                case ImportMode.Boks:
                    AddNumberInput("Høyde", "H", "Boksens høyde (Z).", 2.0);
                    AddNumberInput("Bredde", "B", "Boksens bredde (X).", 1.0);
                    AddNumberInput("Lengde", "L", "Boksens lengde (Y).", 1.0);
                    break;
            }

            Params.OnParametersChanged();
        }

        private void EnsureOutputsForMode(ImportMode mode, List<string> typeKeysForGenerell)
        {
            // Clear all outputs
            while (Params.Output.Count > 0)
                Params.UnregisterOutputParameter(Params.Output[0], true);

            if (mode == ImportMode.Generell)
            {
                AddPointOutput("Punkt", "Pts", "Alle .PUNKT (raw).");
                AddCurveOutput("Kurver", "Crv", "Alle .KURVE som polyline (raw).");
                AddBrepOutput("Flater", "Srf", "Alle .FLATE som planar Brep hvis mulig.");

                // OBJTYPE buckets (Param_Geometry)
                var types = (typeKeysForGenerell == null || typeKeysForGenerell.Count == 0)
                    ? new List<string> { "Udefinert" }
                    : typeKeysForGenerell;

                var usedNicks = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

                foreach (var t in types)
                {
                    var name = $"OBJTYPE: {t}";
                    var nick = MakeUniqueNick(SanitizeNick(t), usedNicks);
                    AddGeometryOutput(name, nick, $"Alle objekter (punkt/kurve/flate) med OBJTYPE '{t}'.");
                }

                AddTextOutput("Info", "i", "Status og diagnostikk.");
                Params.OnParametersChanged();
                return;
            }

            if (mode == ImportMode.Ror)
            {
                AddBrepOutput("Rør", "Pipes", "Rør-breps laget fra alle .KURVE.");
                AddTextOutput("ObjType", "types", "OBJTYPE per rør (1-til-1). Mangler -> 'Udefinert'.");
                AddTextOutput("Info", "i", "Status og diagnostikk.");
                Params.OnParametersChanged();
                return;
            }

            if (mode == ImportMode.Boks)
            {
                AddBrepOutput("Bokser", "Boxes", "Bokser laget fra alle .PUNKT.");
                AddTextOutput("ObjType", "types", "OBJTYPE per boks (1-til-1). Mangler -> 'Udefinert'.");
                AddTextOutput("Info", "i", "Status og diagnostikk.");
                Params.OnParametersChanged();
                return;
            }
        }

        private void AddNumberInput(string name, string nick, string desc, double defaultValue)
        {
            var p = new Param_Number
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.item
            };
            p.SetPersistentData(new GH_Number(defaultValue));
            Params.RegisterInputParam(p);
        }

        private void AddPointOutput(string name, string nick, string desc)
        {
            var p = new Param_Point
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.list
            };
            Params.RegisterOutputParam(p);
        }

        private void AddCurveOutput(string name, string nick, string desc)
        {
            var p = new Param_Curve
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.list
            };
            Params.RegisterOutputParam(p);
        }

        private void AddBrepOutput(string name, string nick, string desc)
        {
            var p = new Param_Brep
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.list
            };
            Params.RegisterOutputParam(p);
        }

        private void AddTextOutput(string name, string nick, string desc)
        {
            var p = new Param_String
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.list
            };
            Params.RegisterOutputParam(p);
        }

        private void AddGeometryOutput(string name, string nick, string desc)
        {
            var p = new Param_Geometry
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.list
            };
            Params.RegisterOutputParam(p);
        }

        // IGH_VariableParameterComponent (vi styrer parametre selv)
        public bool CanInsertParameter(GH_ParameterSide side, int index) => false;
        public bool CanRemoveParameter(GH_ParameterSide side, int index) => false;
        public IGH_Param CreateParameter(GH_ParameterSide side, int index) => null;
        public bool DestroyParameter(GH_ParameterSide side, int index) => false;
        public void VariableParameterMaintenance() { }

        // ------------------------ TYPE / SIGNATURE HELPERS ------------------------

        private static string NormalizeType(string objType)
        {
            var t = (objType ?? "").Trim();
            return t.Length == 0 ? "Udefinert" : t;
        }

        private static List<string> CollectAllObjTypes(SosiData d)
        {
            var set = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

            bool anyUndef = false;

            foreach (var p in d.Points.Where(x => x.Pos.IsValid))
            {
                var t = (p.ObjType ?? "").Trim();
                if (t.Length == 0) anyUndef = true;
                else set.Add(t);
            }

            foreach (var c in d.Curves.Where(x => x.Points.Count >= 2))
            {
                var t = (c.ObjType ?? "").Trim();
                if (t.Length == 0) anyUndef = true;
                else set.Add(t);
            }

            foreach (var f in d.Flater.Where(x => x.Points.Count >= 3))
            {
                var t = (f.ObjType ?? "").Trim();
                if (t.Length == 0) anyUndef = true;
                else set.Add(t);
            }

            // Hvis ingen typer i det hele tatt: vi vil fortsatt ha Udefinert-bucket
            if (set.Count == 0) anyUndef = true;

            var list = set.OrderBy(s => s, StringComparer.OrdinalIgnoreCase).ToList();

            // legg Udefinert først (kan endres hvis du vil)
            if (anyUndef)
                list.Insert(0, "Udefinert");

            return list;
        }

        private static string SanitizeNick(string s)
        {
            var t = NormalizeType(s);
            var chars = t.Select(ch => char.IsLetterOrDigit(ch) ? ch : '_').ToArray();
            var nick = new string(chars);

            if (string.IsNullOrWhiteSpace(nick))
                nick = "Udefinert";

            if (char.IsDigit(nick[0]))
                nick = "T_" + nick;

            // begrens litt (GH liker korte nicknames)
            if (nick.Length > 24)
                nick = nick.Substring(0, 24);

            return nick;
        }

        private static string MakeUniqueNick(string baseNick, HashSet<string> used)
        {
            var nick = baseNick;
            int i = 2;
            while (used.Contains(nick))
            {
                var suffix = "_" + i;
                var cut = Math.Max(1, 24 - suffix.Length);
                nick = (baseNick.Length > cut ? baseNick.Substring(0, cut) : baseNick) + suffix;
                i++;
            }
            used.Add(nick);
            return nick;
        }

        // ------------------------ FLATE -> BREP ------------------------

        private static Brep TryCreatePlanarBrep(List<Point3d> pts)
        {
            try
            {
                if (pts == null || pts.Count < 3) return null;

                var pl = new Polyline(pts);
                if (!pl.IsClosed) pl.Add(pl[0]);

                var crv = new PolylineCurve(pl);
                if (!crv.IsClosed) return null;

                var tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
                var breps = Brep.CreatePlanarBreps(crv, tol);
                return (breps != null && breps.Length > 0) ? breps[0] : null;
            }
            catch
            {
                return null;
            }
        }

        // ------------------------ PARSER ------------------------

        private sealed class SosiData
        {
            public double Enhet { get; set; } = 1.0;
            public double OrigoN { get; set; } = 0.0;
            public double OrigoO { get; set; } = 0.0;

            public List<SosiPoint> Points { get; } = new List<SosiPoint>();
            public List<SosiCurve> Curves { get; } = new List<SosiCurve>();
            public List<SosiFlate> Flater { get; } = new List<SosiFlate>();
        }

        private sealed class SosiPoint
        {
            public string ObjType { get; set; } = "";
            public Point3d Pos { get; set; } = Point3d.Unset;
        }

        private sealed class SosiCurve
        {
            public string ObjType { get; set; } = "";
            public List<Point3d> Points { get; } = new List<Point3d>();
        }

        private sealed class SosiFlate
        {
            public string ObjType { get; set; } = "";
            public List<Point3d> Points { get; } = new List<Point3d>();
        }

        private static class SosiParser
        {
            public static SosiData Parse(string[] lines)
            {
                var d = new SosiData();
                var inv = CultureInfo.InvariantCulture;

                bool inHode = false;

                SosiPoint currentPoint = null;
                SosiCurve currentCurve = null;
                SosiFlate currentFlate = null;

                foreach (var raw in lines)
                {
                    var line = raw.Trim();
                    if (line.Length == 0) continue;

                    // Block starts
                    if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                    {
                        inHode = true;
                        currentPoint = null;
                        currentCurve = null;
                        currentFlate = null;
                        continue;
                    }

                    if (line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                    {
                        inHode = false;
                        currentPoint = new SosiPoint();
                        currentCurve = null;
                        currentFlate = null;
                        d.Points.Add(currentPoint);
                        continue;
                    }

                    if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                    {
                        inHode = false;
                        currentCurve = new SosiCurve();
                        currentPoint = null;
                        currentFlate = null;
                        d.Curves.Add(currentCurve);
                        continue;
                    }

                    if (line.StartsWith(".FLATE", StringComparison.OrdinalIgnoreCase))
                    {
                        inHode = false;
                        currentFlate = new SosiFlate();
                        currentPoint = null;
                        currentCurve = null;
                        d.Flater.Add(currentFlate);
                        continue;
                    }

                    // Header parsing
                    if (inHode)
                    {
                        if (line.IndexOf("ORIGO", StringComparison.OrdinalIgnoreCase) >= 0)
                        {
                            var nums = SplitNumbers(line);
                            if (nums.Count >= 2)
                            {
                                if (TryParseDouble(nums[0], out var n)) d.OrigoN = n;
                                if (TryParseDouble(nums[1], out var o)) d.OrigoO = o;
                            }
                            continue;
                        }

                        if (line.IndexOf("ENHET", StringComparison.OrdinalIgnoreCase) >= 0)
                        {
                            var nums = SplitNumbers(line);
                            if (nums.Count >= 1 && TryParseDouble(nums[0], out var e)) d.Enhet = e;
                            continue;
                        }

                        continue;
                    }

                    // OBJTYPE
                    if (line.StartsWith("..OBJTYPE", StringComparison.OrdinalIgnoreCase))
                    {
                        var val = line.Substring("..OBJTYPE".Length).Trim().Trim(':').Trim();
                        if (currentPoint != null) currentPoint.ObjType = val;
                        if (currentCurve != null) currentCurve.ObjType = val;
                        if (currentFlate != null) currentFlate.ObjType = val;
                        continue;
                    }

                    // Skip some non-coordinate markers
                    if (line.StartsWith("..NØH", StringComparison.OrdinalIgnoreCase))
                        continue;

                    // Coordinate line: N O H (N=Nord, O=Øst, H=høyde)
                    var nums2 = SplitNumbers(line);
                    if (nums2.Count >= 3 &&
                        double.TryParse(nums2[0], NumberStyles.Float, inv, out double N) &&
                        double.TryParse(nums2[1], NumberStyles.Float, inv, out double O) &&
                        double.TryParse(nums2[2], NumberStyles.Float, inv, out double H))
                    {
                        var pt = new Point3d(
                            d.OrigoO + O * d.Enhet,
                            d.OrigoN + N * d.Enhet,
                            H * d.Enhet);

                        if (currentPoint != null)
                        {
                            // Punkt har typisk én koordinatlinje
                            if (!currentPoint.Pos.IsValid)
                                currentPoint.Pos = pt;
                        }
                        else if (currentCurve != null)
                        {
                            currentCurve.Points.Add(pt);
                        }
                        else if (currentFlate != null)
                        {
                            currentFlate.Points.Add(pt);
                        }
                    }
                }

                // Fjern punkter uten gyldig posisjon
                d.Points.RemoveAll(p => !p.Pos.IsValid);

                return d;
            }

            private static bool TryParseDouble(string s, out double v)
            {
                return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v)
                    || double.TryParse(s.Replace(',', '.'), NumberStyles.Float, CultureInfo.InvariantCulture, out v);
            }

            private static List<string> SplitNumbers(string line)
            {
                var res = new List<string>();
                if (string.IsNullOrWhiteSpace(line)) return res;

                var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                foreach (var p in parts)
                {
                    var s = p.Replace(',', '.');
                    if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
                        res.Add(s);
                }
                return res;
            }
        }
    }
}