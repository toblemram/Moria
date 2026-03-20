using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;

namespace Moria.ReFac
{
    /// <summary>
    /// SOSI.Import (Én komponent)
    ///
    /// Modus:
    /// 0 = Generell
    ///   - Outputs: Punkt (alle .PUNKT), Kurver (alle .KURVE), Flater (alle .FLATE som planar brep hvis mulig)
    ///   - + egne outputs per OBJTYPE som inneholder alle objekter (punkt/kurve/flate)
    ///   - Manglende OBJTYPE -> "Udefinert"
    ///
    /// 1 = Rør
    ///   - Input: Diameter
    ///   - Outputs: Pipes (Brep), Types (tekst 1-til-1), Info
    ///
    /// 2 = Boks
    ///   - Inputs: Høyde, Bredde, Lengde, (valgfri) Senterlinje
    ///   - Outputs: Bokser (Brep), Types (tekst 1-til-1), Info
    ///   - Hvis senterlinje er satt: boksene roteres rundt Z etter tangent ved nærmeste punkt (tangent i XY).
    ///
    /// 3 = Skilt
    ///   - Inputs: Radius, (valgfri) Senterlinje
    ///   - Leser PUNKT. I hvert punkt lages en oppreist sirkel (Curve) orientert vinkelrett på senterlinje.
    ///   - Outputs: Skilt (Curve), Types (tekst 1-til-1), Info.
    ///
    /// 4 = Vifte
    ///   - Inputs: Radius, Avstand, (valgfri) Senterlinje
    ///   - Leser PUNKT. I hvert punkt lages en "oppreist" sirkel der punktet er bunnen (senter = punkt + Z*radius).
    ///   - Sirkelen orienteres vinkelrett på senterlinjen: plan normal ~ tangent (projisert i XY).
    ///   - Sirkelen kopieres og flyttes Avstand frem og tilbake langs tangent-retning, og loftes til en sylinder (Brep).
    ///   - Outputs: Sylindre (Brep), Types (tekst 1-til-1), Info.
    /// </summary>
    public class GH_SosiImport : GH_Component, IGH_VariableParameterComponent
    {
        public GH_SosiImport() : base(
            "SOSI.Import",
            "SOSI",
            "Leser SOSI. Modus: Generell, Rør, Boks, Skilt, Vifte.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("2E2C3D7E-5C8E-4A19-9D3F-0E6F1C8E7C12");
        protected override System.Drawing.Bitmap Icon => null;

        private enum Mode { Generell = 0, Ror = 1, Boks = 2, Skilt = 3, Vifte = 4 }

        private Mode _lastMode = (Mode)(-1);
        private string _lastTypeSig = "";

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);

            p.AddIntegerParameter("Mode", "mode", "0=Generell, 1=Rør, 2=Boks, 3=Skilt, 4=Vifte", GH_ParamAccess.item, 0);
            if (p[1] is Param_Integer m)
            {
                m.AddNamedValue("Generell", (int)Mode.Generell);
                m.AddNamedValue("Rør", (int)Mode.Ror);
                m.AddNamedValue("Boks", (int)Mode.Boks);
                m.AddNamedValue("Skilt", (int)Mode.Skilt);
                m.AddNamedValue("Vifte", (int)Mode.Vifte);
            }

            EnsureInputsForMode(Mode.Generell);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            // Start-sett (Generell). OBJTYPE-buckets bygges dynamisk etter parsing.
            p.AddPointParameter("Punkt", "Pts", "Alle .PUNKT", GH_ParamAccess.list);
            p.AddCurveParameter("Kurver", "Crv", "Alle .KURVE", GH_ParamAccess.list);
            p.AddBrepParameter("Flater", "Srf", "Alle .FLATE (planar brep hvis mulig)", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Status", GH_ParamAccess.list);
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

            var mode = (Mode)Math.Max(0, Math.Min(4, modeInt));

            // Mode endret -> rebuild inputs/outputs og re-solve
            if (mode != _lastMode)
            {
                _lastMode = mode;
                EnsureInputsForMode(mode);
                EnsureOutputsForMode(mode, typeKeysForGenerell: new List<string> { "Udefinert" });
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

            if (mode == Mode.Generell)
            {
                var typeKeys = SosiUtils.CollectAllObjTypes(data);
                var sig = string.Join("|", typeKeys).ToLowerInvariant();

                if (!string.Equals(sig, _lastTypeSig, StringComparison.Ordinal))
                {
                    _lastTypeSig = sig;
                    EnsureOutputsForMode(mode, typeKeys);
                    Params.OnParametersChanged();
                    ExpireSolution(true);
                    return;
                }

                SolveGenerell(da, data, typeKeys);
                return;
            }

            if (mode == Mode.Ror)
            {
                double diameter = 1.0;
                da.GetData(2, ref diameter);

                var result = SosiGeometry.BuildPipes(data, diameter);
                da.SetDataList(0, result.Breps);
                da.SetDataList(1, result.Types);
                da.SetDataList(2, result.Info);
                return;
            }

            if (mode == Mode.Boks)
            {
                double h = 2.0, b = 1.0, l = 1.0;
                da.GetData(2, ref h);
                da.GetData(3, ref b);
                da.GetData(4, ref l);

                Curve centerline = null;
                bool hasCenter = false;
                if (Params.Input.Count >= 6) // siste er optional
                    hasCenter = da.GetData(5, ref centerline);

                var result = SosiGeometry.BuildBoxes(data, h, b, l, hasCenter ? centerline : null);
                da.SetDataList(0, result.Breps);
                da.SetDataList(1, result.Types);
                da.SetDataList(2, result.Info);
                return;
            }

            if (mode == Mode.Skilt)
            {
                double radius = 1.0;
                da.GetData(2, ref radius);

                Curve centerline = null;
                bool hasCenter = false;
                if (Params.Input.Count >= 4) // optional
                    hasCenter = da.GetData(3, ref centerline);

                var result = SosiGeometry.BuildSigns(data, radius, hasCenter ? centerline : null);
                da.SetDataList(0, result.Curves);
                da.SetDataList(1, result.Types);
                da.SetDataList(2, result.Info);
                return;
            }

            // Vifte
            {
                double radius = 1.0;
                double dist = 1.0;
                da.GetData(2, ref radius);
                da.GetData(3, ref dist);

                Curve centerline = null;
                bool hasCenter = false;
                if (Params.Input.Count >= 5) // optional
                    hasCenter = da.GetData(4, ref centerline);

                var result = SosiGeometry.BuildVifteCylinders(data, radius, dist, hasCenter ? centerline : null);
                da.SetDataList(0, result.Breps);
                da.SetDataList(1, result.Types);
                da.SetDataList(2, result.Info);
            }
        }

        private void SolveGenerell(IGH_DataAccess da, SosiData data, List<string> typeKeys)
        {
            var pts = SosiGeometry.MakePointList(data);
            var curves = SosiGeometry.MakeCurveList(data);
            var flater = SosiGeometry.MakeFlateBrepList(data);

            var buckets = SosiGeometry.BuildObjTypeBuckets(
                data,
                flaterByFlateIndex: flater.FlaterByIndex,
                typeKeys: typeKeys);

            var info = new List<string>
            {
                "Mode=Generell",
                $"PUNKT={pts.Count}",
                $"KURVE={curves.Count}",
                $".FLATE(brep)={flater.Breps.Count} (feil={flater.FailCount})",
                $"OBJTYPE buckets={typeKeys.Count}",
                $"ENHET={data.Enhet}, ORIGO(N,O)=({data.OrigoN},{data.OrigoO})"
            };

            da.SetDataList(0, pts);
            da.SetDataList(1, curves);
            da.SetDataList(2, flater.Breps);

            int start = 3;
            for (int i = 0; i < typeKeys.Count; i++)
            {
                var key = typeKeys[i];
                da.SetDataList(start + i, buckets.TryGetValue(key, out var list) ? list : new List<object>());
            }

            da.SetDataList(start + typeKeys.Count, info);
        }

        private void EnsureInputsForMode(Mode mode)
        {
            // behold alltid [0]=path og [1]=mode
            while (Params.Input.Count > 2)
                Params.UnregisterInputParameter(Params.Input[2], true);

            switch (mode)
            {
                case Mode.Generell:
                    break;

                case Mode.Ror:
                    GhParamUtils.AddNumberInput(Params, "Diameter", "D", "Rørdiameter for alle KURVE.", 1.0);
                    break;

                case Mode.Boks:
                    GhParamUtils.AddNumberInput(Params, "Høyde", "H", "Boksens høyde (Z).", 2.0);
                    GhParamUtils.AddNumberInput(Params, "Bredde", "B", "Boksens bredde (X).", 1.0);
                    GhParamUtils.AddNumberInput(Params, "Lengde", "L", "Boksens lengde (Y).", 1.0);
                    GhParamUtils.AddOptionalCurveInput(Params, "Senterlinje", "CL",
                        "Valgfri kurve. Roterer boksene etter nærmeste punkt (tangent i XY).");
                    break;

                case Mode.Skilt:
                    GhParamUtils.AddNumberInput(Params, "Radius", "R", "Radius på skilt-sirkler.", 1.0);
                    GhParamUtils.AddOptionalCurveInput(Params, "Senterlinje", "CL",
                        "Valgfri kurve. Orienterer skilt vinkelrett på senterlinjen (tangent i XY).");
                    break;

                case Mode.Vifte:
                    GhParamUtils.AddNumberInput(Params, "Radius", "R", "Radius på sylinder (basert på sirkel).", 1.0);
                    GhParamUtils.AddNumberInput(Params, "Avstand", "S", "Avstand frem/tilbake langs tangent før loft (total lengde ~ 2*S).", 1.0);
                    GhParamUtils.AddOptionalCurveInput(Params, "Senterlinje", "CL",
                        "Valgfri kurve. Orienterer sirkel vinkelrett på senterlinjen (tangent i XY) og flytter langs tangent.");
                    break;
            }
        }

        private void EnsureOutputsForMode(Mode mode, List<string> typeKeysForGenerell)
        {
            while (Params.Output.Count > 0)
                Params.UnregisterOutputParameter(Params.Output[0], true);

            if (mode == Mode.Generell)
            {
                GhParamUtils.AddPointOutput(Params, "Punkt", "Pts", "Alle .PUNKT");
                GhParamUtils.AddCurveOutput(Params, "Kurver", "Crv", "Alle .KURVE");
                GhParamUtils.AddBrepOutput(Params, "Flater", "Srf", "Alle .FLATE (planar brep hvis mulig)");

                var used = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
                foreach (var t in (typeKeysForGenerell ?? new List<string> { "Udefinert" }))
                {
                    var name = $"OBJTYPE: {t}";
                    var nick = GhParamUtils.MakeUniqueNick(GhParamUtils.SanitizeNick(t), used);
                    GhParamUtils.AddGeometryOutput(Params, name, nick,
                        $"Alle objekter (punkt/kurve/flate) med OBJTYPE '{t}'.");
                }

                GhParamUtils.AddTextOutput(Params, "Info", "i", "Status");
                return;
            }

            if (mode == Mode.Ror)
            {
                GhParamUtils.AddBrepOutput(Params, "Rør", "Pipes", "Pipes fra alle KURVE");
                GhParamUtils.AddTextOutput(Params, "ObjType", "types", "OBJTYPE per rør (1-til-1). Mangler -> Udefinert.");
                GhParamUtils.AddTextOutput(Params, "Info", "i", "Status");
                return;
            }

            if (mode == Mode.Boks)
            {
                GhParamUtils.AddBrepOutput(Params, "Bokser", "Boxes", "Bokser fra alle PUNKT");
                GhParamUtils.AddTextOutput(Params, "ObjType", "types", "OBJTYPE per boks (1-til-1). Mangler -> Udefinert.");
                GhParamUtils.AddTextOutput(Params, "Info", "i", "Status");
                return;
            }

            if (mode == Mode.Skilt)
            {
                GhParamUtils.AddCurveOutput(Params, "Skilt", "Signs", "Oppreiste sirkler (Curve) per PUNKT, vinkelrett på senterlinje.");
                GhParamUtils.AddTextOutput(Params, "ObjType", "types", "OBJTYPE per skilt (1-til-1). Mangler -> Udefinert.");
                GhParamUtils.AddTextOutput(Params, "Info", "i", "Status");
                return;
            }

            // Vifte
            GhParamUtils.AddBrepOutput(Params, "Sylindre", "Cyl", "Sylinder-breps fra PUNKT (loft av to sirkler flyttet ±Avstand).");
            GhParamUtils.AddTextOutput(Params, "ObjType", "types", "OBJTYPE per sylinder (1-til-1). Mangler -> Udefinert.");
            GhParamUtils.AddTextOutput(Params, "Info", "i", "Status");
        }

        // IGH_VariableParameterComponent (vi styrer selv)
        public bool CanInsertParameter(GH_ParameterSide side, int index) => false;
        public bool CanRemoveParameter(GH_ParameterSide side, int index) => false;
        public IGH_Param CreateParameter(GH_ParameterSide side, int index) => null;
        public bool DestroyParameter(GH_ParameterSide side, int index) => false;
        public void VariableParameterMaintenance() { }
    }
}
