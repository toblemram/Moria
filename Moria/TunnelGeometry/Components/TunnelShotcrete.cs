using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;

using Rhino;
using Rhino.Geometry;

using Moria.TunnelGeometry; // ShotcreteProfiles2D

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// Sweeps zero-thickness shotcrete + insulation profiles along the path
    /// and thickens them. Also generates rebar grid as curves + solid pipes.
    ///
    /// Geometry:
    ///  - Shotcrete:
    ///      * 2D open profile langs T-profilen med ett skjørt.
    ///      * Sweeps langs path → flate.
    ///      * Offset til solid med tykkelse tShot.
    ///  - Insulation:
    ///      * Bygges utenpå shotcrete (betongtykkelse mellom), samme form,
    ///        men med ekstra skjørt nr. 2.
    ///  - Rebar:
    ///      * Ringer (tverrsnitt) med spacing sL langs tunnelen.
    ///      * Langsgående stenger rundt tverrsnittet med spacing sT.
    ///      * Rebarposisjon velges fra shotcrete-profilet, med cover ut i betongen.
    ///      * Rebar følger base- / bay-profil avhengig av stasjon (busslommer).
    /// </summary>
    public class GH_TunnelShotcrete3D : GH_Component
    {
        public GH_TunnelShotcrete3D()
            : base(
                "Tunnel Shotcrete 3D",
                "Shotcrete3D",
                "Shotcrete + insulation + rebar along a tunnel path, including emergency bays.",
                "Tunnel",
                "Tunnel Geom")
        {
        }

        // --------------------------------------------------------------------
        // INPUTS
        // --------------------------------------------------------------------

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddCurveParameter(
                "Path",
                "P",
                "Tunnel centreline/path.",
                GH_ParamAccess.item);

            p.AddTextParameter(
                "BaseProfileType",
                "T",
                "Base tunnel T-profile type, e.g. T5.5, T7.5, T8.5, T9.5, T10.5, T12.5, T13, T13.5, T14.",
                GH_ParamAccess.item,
                "T10.5");

            p.AddNumberParameter(
                "ShotcreteThickness",
                "tShot",
                "Shotcrete thickness [m], applied as an offset from the inner shotcrete surface.",
                GH_ParamAccess.item,
                0.15);

            p.AddNumberParameter(
                "SkirtStartY",
                "Yk",
                "Height [m] above the tunnel bottom where the shotcrete skirt starts.",
                GH_ParamAccess.item,
                1.5);

            p.AddNumberParameter(
                "SkirtAngle",
                "Ang",
                "Shotcrete skirt angle [deg] measured from vertical downward direction, bending outwards.",
                GH_ParamAccess.item,
                20.0);

            p.AddNumberParameter(
                "SkirtLength",
                "L",
                "Shotcrete skirt length [m] from SkirtStartY along the given angle.",
                GH_ParamAccess.item,
                2.0);

            // Isolasjon
            p.AddNumberParameter(
                "InsulThickness",
                "tIso",
                "Insulation thickness [m] outside the shotcrete.",
                GH_ParamAccess.item,
                0.10);

            p.AddNumberParameter(
                "InsulSkirtAngle",
                "AngIso",
                "Insulation second skirt angle [deg] from vertical, starting at the toe of the shotcrete outer skirt.",
                GH_ParamAccess.item,
                10.0);

            p.AddNumberParameter(
                "InsulSkirtLength",
                "LIso",
                "Insulation second skirt length [m].",
                GH_ParamAccess.item,
                1.0);

            // Busslommer
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

            // Armering
            p.AddNumberParameter(
                "RebarTransSpacing",
                "sT",
                "Rebar spacing along the cross-section (around the tunnel) [m].",
                GH_ParamAccess.item,
                0.25);

            p.AddNumberParameter(
                "RebarLongSpacing",
                "sL",
                "Rebar spacing along the tunnel [m].",
                GH_ParamAccess.item,
                0.30);

            p.AddNumberParameter(
                "RebarDiameter",
                "Ø",
                "Rebar diameter [m]. Used to create solid pipes.",
                GH_ParamAccess.item,
                0.016);

            p.AddNumberParameter(
                "RebarCover",
                "c",
                "Concrete cover [m] from inner shotcrete surface out into the concrete to the rebar grid.",
                GH_ParamAccess.item,
                0.05);
        }

        // --------------------------------------------------------------------
        // OUTPUTS
        // --------------------------------------------------------------------

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter(
                "Shotcrete",
                "Shot",
                "Shotcrete tunnel Brep (solid, can be used for volume).",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Shotcrete Segment Solids",
                "ShotSeg",
                "Shotcrete solids per segment (debug).",
                GH_ParamAccess.list);

            p.AddBrepParameter(
                "Insulation",
                "Iso",
                "Insulation Brep (solid) outside the shotcrete.",
                GH_ParamAccess.item);

            p.AddBrepParameter(
                "Insulation Segment Solids",
                "IsoSeg",
                "Insulation solids per segment (debug).",
                GH_ParamAccess.list);

            p.AddCurveParameter(
                "Rebar Curves",
                "Rebar",
                "Rebar curves (transverse rings + longitudinal bars).",
                GH_ParamAccess.list);

            p.AddBrepParameter(
                "Rebar Solids",
                "RebarB",
                "Rebar as solid pipes with the specified diameter.",
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

        private class BayDef
        {
            public double SCenter;
            public int Side;      // +1 = right, -1 = left
            public double S0;     // start transition up
            public double S1;     // start full bay
            public double S2;     // end   full bay
            public double S3;     // end   transition down
        }

        private class SegmentDef
        {
            public double SStart;
            public double SEnd;

            public Curve InnerShotStart2D;
            public Curve InnerShotEnd2D;

            public Curve InnerIsoStart2D;
            public Curve InnerIsoEnd2D;
        }

        // --------------------------------------------------------------------
        // MAIN
        // --------------------------------------------------------------------

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Curve path = null;
            string baseProfileType = "T10.5";
            double tShot = 0.15;
            double skirtStartY = 1.5;
            double skirtAngleDeg = 20.0;
            double skirtLength = 2.0;

            double tIso = 0.10;
            double isoAngDeg = 10.0;
            double isoLen = 1.0;

            var stationsRight = new List<double>();
            var stationsLeft = new List<double>();

            double sT = 0.25;
            double sL = 0.30;
            double rebarDia = 0.016;
            double rebarCover = 0.05;

            da.GetData(0, ref path);
            da.GetData(1, ref baseProfileType);
            da.GetData(2, ref tShot);
            da.GetData(3, ref skirtStartY);
            da.GetData(4, ref skirtAngleDeg);
            da.GetData(5, ref skirtLength);
            da.GetData(6, ref tIso);
            da.GetData(7, ref isoAngDeg);
            da.GetData(8, ref isoLen);
            da.GetDataList(9, stationsRight);
            da.GetDataList(10, stationsLeft);
            da.GetData(11, ref sT);
            da.GetData(12, ref sL);
            da.GetData(13, ref rebarDia);
            da.GetData(14, ref rebarCover);

            var info = new List<string>();
            var shotcreteSegments = new List<Brep>();
            var insulationSegments = new List<Brep>();
            var rebarCurves = new List<Curve>();
            var rebarSolids = new List<Brep>();

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

            if (tShot <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ShotcreteThickness must be > 0.");
                return;
            }

            if (skirtLength <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SkirtLength must be > 0.");
                return;
            }

            if (sT <= 0.0 || sL <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Rebar spacings must be > 0.");
                return;
            }

            if (rebarDia <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "RebarDiameter must be > 0.");
                return;
            }

            baseProfileType = baseProfileType.Replace(",", ".").ToUpperInvariant();

            info.Add($"Path length: {pathLength:0.###} m");
            info.Add($"BaseProfileType: {baseProfileType}");
            info.Add($"Shotcrete t: {tShot:0.###} m");
            info.Add($"SkirtStartY: {skirtStartY:0.###} m");
            info.Add($"SkirtAngle: {skirtAngleDeg:0.###}°");
            info.Add($"SkirtLength: {skirtLength:0.###} m");
            info.Add($"Insul t: {tIso:0.###} m (0 = off)");
            info.Add($"Insul skirt angle: {isoAngDeg:0.###}°");
            info.Add($"Insul skirt length: {isoLen:0.###} m");
            info.Add($"Rebar sT: {sT:0.###} m");
            info.Add($"Rebar sL: {sL:0.###} m");
            info.Add($"Rebar Ø: {rebarDia * 1000.0:0.#} mm");
            info.Add($"Rebar cover (outwards): {rebarCover:0.###} m");
            info.Add($"Model tol: {tol:0.###E0}");

            // ---------------------------------------------------------------
            // 1) Profile parameters
            // ---------------------------------------------------------------

            if (!ProfileType.Profiles.TryGetValue(baseProfileType, out var basePar))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Base profile type '{baseProfileType}' not found. Valid keys: {string.Join(", ", ProfileType.Profiles.Keys)}");
                return;
            }

            string bayProfileType = ResolveBayProfileType(baseProfileType);
            if (!ProfileType.Profiles.TryGetValue(bayProfileType, out var bayPar))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Bay profile type '{bayProfileType}' not found.");
                return;
            }

            info.Add($"BayProfileType: {bayProfileType}");

            // ---------------------------------------------------------------
            // 2) Build tunnel profiles (2D)
            // ---------------------------------------------------------------

            if (!BuildTunnelProfileWorld(
                    baseProfileType, basePar, tol,
                    out Curve baseTunnelWorld,
                    out string errBase))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBase);
                return;
            }

            if (!BuildTunnelProfileWorld(
                    bayProfileType, bayPar, tol,
                    out Curve bayTunnelWorld,
                    out string errBay))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, errBay);
                return;
            }

            BoundingBox bbBase = baseTunnelWorld.GetBoundingBox(true);
            BoundingBox bbBay = bayTunnelWorld.GetBoundingBox(true);

            double xMinBase = bbBase.Min.X;
            double xMaxBase = bbBase.Max.X;
            double xMinBay = bbBay.Min.X;
            double xMaxBay = bbBay.Max.X;

            double dxRight = xMinBase - xMinBay;
            double dxLeft = xMaxBase - xMaxBay;

            Curve bayTunnelRightWorld = bayTunnelWorld.DuplicateCurve();
            bayTunnelRightWorld.Transform(Transform.Translation(dxRight, 0, 0));

            Curve bayTunnelLeftWorld = bayTunnelWorld.DuplicateCurve();
            bayTunnelLeftWorld.Transform(Transform.Translation(dxLeft, 0, 0));

            info.Add($"Base X-span: [{xMinBase:0.###}, {xMaxBase:0.###}]");
            info.Add($"Bay  X-span: [{xMinBay:0.###}, {xMaxBay:0.###}]");
            info.Add($"Right bay dx = {dxRight:0.###}");
            info.Add($"Left  bay dx = {dxLeft:0.###}");

            // ---------------------------------------------------------------
            // 3) Build shotcrete inner profiles (2D net-profiler)
            // ---------------------------------------------------------------

            if (!ShotcreteProfiles2D.BuildShotcreteInnerProfile2D(
                    baseTunnelWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve baseShotInner2D,
                    out string scErrBase))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Base shotcrete: " + scErrBase);
                return;
            }

            if (!ShotcreteProfiles2D.BuildShotcreteInnerProfile2D(
                    bayTunnelRightWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve bayShotInnerRight2D,
                    out string scErrBayR))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Right bay shotcrete: " + scErrBayR);
                return;
            }

            if (!ShotcreteProfiles2D.BuildShotcreteInnerProfile2D(
                    bayTunnelLeftWorld,
                    skirtStartY,
                    skirtAngleDeg,
                    skirtLength,
                    tol,
                    out Curve bayShotInnerLeft2D,
                    out string scErrBayL))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Left bay shotcrete: " + scErrBayL);
                return;
            }

            // ---------------------------------------------------------------
            // 4) Build insulation inner profiles (2D), based on shotcrete
            // ---------------------------------------------------------------

            Curve baseIsoInner2D = null;
            Curve bayIsoInnerRight2D = null;
            Curve bayIsoInnerLeft2D = null;

            if (tIso > 0.0)
            {
                if (!ShotcreteProfiles2D.BuildInsulationInnerProfile2D(
                        baseShotInner2D,
                        tShot,
                        isoAngDeg,
                        isoLen,
                        tol,
                        out baseIsoInner2D,
                        out string isoErrBase))
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Base insulation: " + isoErrBase);
                    tIso = 0.0;
                }
                else if (!ShotcreteProfiles2D.BuildInsulationInnerProfile2D(
                        bayShotInnerRight2D,
                        tShot,
                        isoAngDeg,
                        isoLen,
                        tol,
                        out bayIsoInnerRight2D,
                        out string isoErrBayR))
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Right bay insulation: " + isoErrBayR);
                    tIso = 0.0;
                }
                else if (!ShotcreteProfiles2D.BuildInsulationInnerProfile2D(
                        bayShotInnerLeft2D,
                        tShot,
                        isoAngDeg,
                        isoLen,
                        tol,
                        out bayIsoInnerLeft2D,
                        out string isoErrBayL))
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Left bay insulation: " + isoErrBayL);
                    tIso = 0.0;
                }
            }

            // ---------------------------------------------------------------
            // 5) Bay definitions
            // ---------------------------------------------------------------

            const double rampLen = 30.0;
            const double flatLen = 30.0;

            var bays = new List<BayDef>();

            foreach (double s in stationsRight)
            {
                bays.Add(new BayDef
                {
                    SCenter = s,
                    Side = +1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                });
            }

            foreach (double s in stationsLeft)
            {
                bays.Add(new BayDef
                {
                    SCenter = s,
                    Side = -1,
                    S0 = s - (rampLen + flatLen / 2.0),
                    S1 = s - (flatLen / 2.0),
                    S2 = s + (flatLen / 2.0),
                    S3 = s + (rampLen + flatLen / 2.0)
                });
            }

            bays.Sort((a, b) => a.SCenter.CompareTo(b.SCenter));

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

            info.Add($"Number of bays: {bays.Count}");

            // ---------------------------------------------------------------
            // 6) Build sweep segments (shotcrete + insulation)
            // ---------------------------------------------------------------

            var segments = new List<SegmentDef>();
            double currentS = 0.0;
            double eps = 1e-6;

            foreach (var b in bays)
            {
                if (b.S0 > currentS + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = currentS,
                        SEnd = b.S0,
                        InnerShotStart2D = baseShotInner2D,
                        InnerShotEnd2D = baseShotInner2D,
                        InnerIsoStart2D = baseIsoInner2D,
                        InnerIsoEnd2D = baseIsoInner2D
                    });
                }

                Curve shotBay2D = (b.Side > 0) ? bayShotInnerRight2D : bayShotInnerLeft2D;
                Curve isoBay2D = null;
                if (tIso > 0.0)
                    isoBay2D = (b.Side > 0) ? bayIsoInnerRight2D : bayIsoInnerLeft2D;

                if (b.S1 > b.S0 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S0,
                        SEnd = b.S1,
                        InnerShotStart2D = baseShotInner2D,
                        InnerShotEnd2D = shotBay2D,
                        InnerIsoStart2D = baseIsoInner2D,
                        InnerIsoEnd2D = isoBay2D
                    });
                }

                if (b.S2 > b.S1 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S1,
                        SEnd = b.S2,
                        InnerShotStart2D = shotBay2D,
                        InnerShotEnd2D = shotBay2D,
                        InnerIsoStart2D = isoBay2D,
                        InnerIsoEnd2D = isoBay2D
                    });
                }

                if (b.S3 > b.S2 + eps)
                {
                    segments.Add(new SegmentDef
                    {
                        SStart = b.S2,
                        SEnd = b.S3,
                        InnerShotStart2D = shotBay2D,
                        InnerShotEnd2D = baseShotInner2D,
                        InnerIsoStart2D = isoBay2D,
                        InnerIsoEnd2D = baseIsoInner2D
                    });
                }

                currentS = Math.Max(currentS, b.S3);
            }

            if (currentS < pathLength - eps)
            {
                segments.Add(new SegmentDef
                {
                    SStart = currentS,
                    SEnd = pathLength,
                    InnerShotStart2D = baseShotInner2D,
                    InnerShotEnd2D = baseShotInner2D,
                    InnerIsoStart2D = baseIsoInner2D,
                    InnerIsoEnd2D = baseIsoInner2D
                });
            }

            if (bays.Count == 0 && segments.Count == 0)
            {
                segments.Add(new SegmentDef
                {
                    SStart = 0.0,
                    SEnd = pathLength,
                    InnerShotStart2D = baseShotInner2D,
                    InnerShotEnd2D = baseShotInner2D,
                    InnerIsoStart2D = baseIsoInner2D,
                    InnerIsoEnd2D = baseIsoInner2D
                });
            }

            info.Add($"Number of sweep segments: {segments.Count}");

            // ---------------------------------------------------------------
            // 7) Sweep per segment (shotcrete + insulation)
            // ---------------------------------------------------------------

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

                // --- Shotcrete ---
                Curve c0Shot = seg.InnerShotStart2D.DuplicateCurve();
                Curve c1Shot = seg.InnerShotEnd2D.DuplicateCurve();

                Transform toStart = Transform.PlaneToPlane(Plane.WorldXY, frame0);
                Transform toEnd = Transform.PlaneToPlane(Plane.WorldXY, frame1);

                c0Shot.Transform(toStart);
                c1Shot.Transform(toEnd);

                Brep[] innerShotBreps = sweep.PerformSweep(rail, new Curve[] { c0Shot, c1Shot });
                if (innerShotBreps == null || innerShotBreps.Length == 0)
                {
                    info.Add($"Segment {segIndex}: shotcrete SweepOneRail failed.");
                }
                else
                {
                    Brep innerSurf = innerShotBreps[0];
                    try
                    {
                        Brep[] blends, walls;
                        Brep[] offset = Brep.CreateOffsetBrep(
                            innerSurf,
                            tShot,
                            true,
                            true,
                            tol,
                            out blends,
                            out walls);

                        if (offset != null && offset.Length > 0)
                        {
                            Brep solid = PrepareSolid(offset[0], $"ShotSeg {segIndex}", tol, info);
                            shotcreteSegments.Add(solid);
                            info.Add($"Segment {segIndex}: shotcrete sweep + offset solid OK.");
                        }
                        else
                        {
                            info.Add($"Segment {segIndex}: shotcrete CreateOffsetBrep failed; using inner surface only.");
                            shotcreteSegments.Add(innerSurf);
                        }
                    }
                    catch (Exception ex)
                    {
                        info.Add($"Segment {segIndex}: shotcrete offset exception: {ex.Message}");
                        shotcreteSegments.Add(innerSurf);
                    }
                }

                // --- Insulation ---
                if (tIso > 0.0 && seg.InnerIsoStart2D != null && seg.InnerIsoEnd2D != null)
                {
                    Curve c0Iso = seg.InnerIsoStart2D.DuplicateCurve();
                    Curve c1Iso = seg.InnerIsoEnd2D.DuplicateCurve();

                    c0Iso.Transform(toStart);
                    c1Iso.Transform(toEnd);

                    Brep[] innerIsoBreps = sweep.PerformSweep(rail, new Curve[] { c0Iso, c1Iso });
                    if (innerIsoBreps == null || innerIsoBreps.Length == 0)
                    {
                        info.Add($"Segment {segIndex}: insulation SweepOneRail failed.");
                    }
                    else
                    {
                        Brep innerIsoSurf = innerIsoBreps[0];
                        try
                        {
                            Brep[] blendsIso, wallsIso;
                            Brep[] offsetIso = Brep.CreateOffsetBrep(
                                innerIsoSurf,
                                tIso,
                                true,
                                true,
                                tol,
                                out blendsIso,
                                out wallsIso);

                            if (offsetIso != null && offsetIso.Length > 0)
                            {
                                Brep solidIso = PrepareSolid(offsetIso[0], $"IsoSeg {segIndex}", tol, info);
                                insulationSegments.Add(solidIso);
                                info.Add($"Segment {segIndex}: insulation sweep + offset solid OK.");
                            }
                            else
                            {
                                info.Add($"Segment {segIndex}: insulation CreateOffsetBrep failed; using inner surface only.");
                                insulationSegments.Add(innerIsoSurf);
                            }
                        }
                        catch (Exception ex)
                        {
                            info.Add($"Segment {segIndex}: insulation offset exception: {ex.Message}");
                            insulationSegments.Add(innerIsoSurf);
                        }
                    }
                }

                segIndex++;
            }

            if (shotcreteSegments.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No shotcrete segment solids were generated.");
                return;
            }

            // ---------------------------------------------------------------
            // 8) Join segments (shotcrete + insulation)
            // ---------------------------------------------------------------

            Brep shotcreteTunnel;
            if (shotcreteSegments.Count == 1)
            {
                shotcreteTunnel = shotcreteSegments[0];
            }
            else
            {
                Brep[] joined = Brep.JoinBreps(shotcreteSegments, tol);
                shotcreteTunnel = (joined != null && joined.Length > 0) ? joined[0] : shotcreteSegments[0];
            }
            shotcreteTunnel = PrepareSolid(shotcreteTunnel, "ShotcreteTunnel", tol, info);

            Brep insulationTunnel = null;
            if (insulationSegments.Count > 0)
            {
                if (insulationSegments.Count == 1)
                    insulationTunnel = insulationSegments[0];
                else
                {
                    Brep[] joinedIso = Brep.JoinBreps(insulationSegments, tol);
                    insulationTunnel = (joinedIso != null && joinedIso.Length > 0) ? joinedIso[0] : insulationSegments[0];
                }
                insulationTunnel = PrepareSolid(insulationTunnel, "InsulationTunnel", tol, info);
            }

            // ---------------------------------------------------------------
            // 9) Rebar generation (curves + pipes) – følger bays
            // ---------------------------------------------------------------

            GenerateRebarGrid(
        path,
        pathLength,
        segments,
        baseShotInner2D,
        sT,
        sL,
        rebarDia,
        rebarCover,
        tol,
        ref rebarCurves,
        ref rebarSolids,
        info);


            // ---------------------------------------------------------------
            // 10) Outputs
            // ---------------------------------------------------------------

            da.SetData(0, shotcreteTunnel);
            da.SetDataList(1, shotcreteSegments);
            da.SetData(2, insulationTunnel);
            da.SetDataList(3, insulationSegments);
            da.SetDataList(4, rebarCurves);
            da.SetDataList(5, rebarSolids);
            da.SetDataList(6, info);
        }

        // --------------------------------------------------------------------
        // REBAR GRID GENERATION
        // --------------------------------------------------------------------
        private static void GenerateRebarGrid(
            Curve path,
            double pathLength,
            List<SegmentDef> segments,
            Curve baseShotInner2D,
            double spacingTransverse,
            double spacingLongitudinal,
            double rebarDia,
            double rebarCover,
            double tol,
            ref List<Curve> rebarCurves,
            ref List<Brep> rebarSolids,
            List<string> info)
        {
            double L = pathLength;

            // ------------------------------------------------------------------
            // Hjelper: finn 2D shotcrete-profil ved stasjon s,
            // ved å bruke segmentene (samme logikk som sweeps).
            // I overgangssegmenter blander vi start- og sluttprofilen lineært.
            // ------------------------------------------------------------------
            Curve GetInnerProfileAtStation(double s)
            {
                foreach (var seg in segments)
                {
                    if (s < seg.SStart - tol)
                        continue;
                    if (s <= seg.SEnd + tol)
                    {
                        Curve a = seg.InnerShotStart2D;
                        Curve b = seg.InnerShotEnd2D;

                        if (a == null && b == null)
                            return baseShotInner2D;
                        if (b == null || seg.SEnd - seg.SStart <= tol || ReferenceEquals(a, b))
                            return a ?? b ?? baseShotInner2D;

                        double u = (s - seg.SStart) / (seg.SEnd - seg.SStart);
                        if (u < 0.0) u = 0.0;
                        if (u > 1.0) u = 1.0;

                        return BlendProfiles2D(a, b, u, tol);
                    }
                }

                return baseShotInner2D;
            }

            // ------------------------------------------------------------------
            // Hjelper: blend mellom to 2D-profiler a og b (i WorldXY) med faktor u.
            // Vi sampler punkter langs begge kurver og interpolerer mellom dem.
            // ------------------------------------------------------------------
            Curve BlendProfiles2D(Curve a, Curve b, double u, double localTol)
            {
                if (a == null && b == null)
                    return null;
                if (a == null) return b.DuplicateCurve();
                if (b == null) return a.DuplicateCurve();

                double lenA = a.GetLength();
                double lenB = b.GetLength();

                if (lenA < localTol || lenB < localTol)
                    return a.DuplicateCurve();

                int n = 80; // antall sample-punkter rundt profilen
                var pts = new List<Point3d>(n + 1);

                for (int i = 0; i <= n; i++)
                {
                    double frac = (double)i / n;

                    double sA = frac * lenA;
                    double sB = frac * lenB;

                    if (!a.LengthParameter(sA, out double tA))
                        continue;
                    if (!b.LengthParameter(sB, out double tB))
                        continue;

                    Point3d pA = a.PointAt(tA);
                    Point3d pB = b.PointAt(tB);

                    Point3d p = (1.0 - u) * pA + u * pB;
                    pts.Add(p);
                }

                if (pts.Count < 2)
                    return a.DuplicateCurve();

                // Glatt kurve gjennom blend-punktene
                return Curve.CreateInterpolatedCurve(pts, 3);
            }

            // ------------------------------------------------------------------
            // Hjelper: offsetter profilet ut i betongen (cover > 0).
            // Vi velger offset-linjen med STØRST bounding box (utover).
            // ------------------------------------------------------------------
            Curve BuildRebar2D(Curve shotInner)
            {
                if (shotInner == null)
                    return null;

                Curve baseCurve = shotInner.DuplicateCurve();
                if (rebarCover <= tol)
                    return baseCurve;

                var candidates = new List<Curve>();

                try
                {
                    Curve[] offPlus = shotInner.Offset(
                        Plane.WorldXY,
                        rebarCover,
                        tol,
                        CurveOffsetCornerStyle.Sharp);
                    if (offPlus != null) candidates.AddRange(offPlus);

                    Curve[] offMinus = shotInner.Offset(
                        Plane.WorldXY,
                        -rebarCover,
                        tol,
                        CurveOffsetCornerStyle.Sharp);
                    if (offMinus != null) candidates.AddRange(offMinus);
                }
                catch
                {
                    return baseCurve;
                }

                Curve best = null;
                double bestDiag = -1.0;
                foreach (var c in candidates)
                {
                    if (c == null) continue;
                    BoundingBox bb = c.GetBoundingBox(true);
                    double diag = bb.Diagonal.Length;
                    if (diag > bestDiag)
                    {
                        bestDiag = diag;
                        best = c;
                    }
                }

                return best ?? baseCurve;
            }

            // ------------------------------------------------------------------
            // 1) Tverrgående ringer (bruker segment-blend for å følge busslommer
            //    OG overgangssoner riktig)
            // ------------------------------------------------------------------

            int nRings = (int)Math.Floor(L / spacingLongitudinal) + 1;
            for (int i = 0; i <= nRings; i++)
            {
                double s = Math.Min(i * spacingLongitudinal, L);

                if (!path.LengthParameter(s, out double t))
                    continue;
                if (!path.PerpendicularFrameAt(t, out Plane frame))
                    continue;

                Curve shotInner = GetInnerProfileAtStation(s);
                if (shotInner == null) continue;

                Curve rebar2D = BuildRebar2D(shotInner);
                if (rebar2D == null) continue;

                Curve ring = rebar2D.DuplicateCurve();
                Transform to3D = Transform.PlaneToPlane(Plane.WorldXY, frame);
                ring.Transform(to3D);
                rebarCurves.Add(ring);

                try
                {
                    double radius = 0.5 * rebarDia;
                    Brep[] pipes = Brep.CreatePipe(
                        ring,
                        radius,
                        false,
                        PipeCapMode.Round,
                        true,
                        tol,
                        tol);
                    if (pipes != null && pipes.Length > 0)
                        rebarSolids.AddRange(pipes);
                }
                catch { }
            }

            // ------------------------------------------------------------------
            // 2) Langsgående stenger (de var allerede "riktige nok", men nå
            //    bruker vi også segment-blend slik at de følger overgangene pent)
            // ------------------------------------------------------------------

            // Referanseprofil ved start (brukes kun for å bestemme hvor mange stenger rundt)
            Curve baseRefShot = GetInnerProfileAtStation(0.0);
            Curve baseRef = BuildRebar2D(baseRefShot);
            if (baseRef == null)
            {
                info?.Add("Rebar: base reference profile is null.");
                return;
            }

            double lenRef = baseRef.GetLength();
            if (lenRef <= tol)
            {
                info?.Add("Rebar: base reference profile length too small.");
                return;
            }

            int nAround = (int)Math.Floor(lenRef / spacingTransverse);
            if (nAround < 1) nAround = 1;

            int samplesAlong = 20;

            for (int idx = 0; idx <= nAround; idx++)
            {
                double uAround = (double)idx / nAround; // 0..1 rundt tverrsnittet
                var pts3D = new List<Point3d>();

                for (int i = 0; i <= samplesAlong; i++)
                {
                    double s = L * (double)i / samplesAlong;
                    if (!path.LengthParameter(s, out double t))
                        continue;
                    if (!path.PerpendicularFrameAt(t, out Plane frame))
                        continue;

                    Curve shotInner = GetInnerProfileAtStation(s);
                    if (shotInner == null) continue;

                    Curve rebar2D = BuildRebar2D(shotInner);
                    if (rebar2D == null) continue;

                    double len = rebar2D.GetLength();
                    if (len < tol) continue;

                    double targetLen = uAround * len;
                    if (!rebar2D.LengthParameter(targetLen, out double tLocal))
                        continue;

                    Point3d p2 = rebar2D.PointAt(tLocal);

                    Point3d p3 = frame.Origin
                                  + frame.XAxis * p2.X
                                  + frame.YAxis * p2.Y;

                    pts3D.Add(p3);
                }

                if (pts3D.Count >= 2)
                {
                    Curve bar = Curve.CreateInterpolatedCurve(pts3D, 3);
                    rebarCurves.Add(bar);

                    try
                    {
                        double radius = 0.5 * rebarDia;
                        Brep[] pipes = Brep.CreatePipe(
                            bar,
                            radius,
                            false,
                            PipeCapMode.Round,
                            true,
                            tol,
                            tol);
                        if (pipes != null && pipes.Length > 0)
                            rebarSolids.AddRange(pipes);
                    }
                    catch { }
                }
            }

            info?.Add($"Rebar curves generated: {rebarCurves.Count}, solids: {rebarSolids.Count}");
        }


        // --------------------------------------------------------------------
        // HELPER: BUILD TUNNEL PROFILE 2D
        // --------------------------------------------------------------------

        private static bool BuildTunnelProfileWorld(
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
            {
                if (string.IsNullOrWhiteSpace(error))
                    error = "Failed to build tunnel profile.";
                return false;
            }

            if (!poly.IsClosed)
            {
                error = "Tunnel profile is not closed.";
                return false;
            }

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
        // BAY PROFILE TYPE HELPER
        // --------------------------------------------------------------------

        private static string ResolveBayProfileType(string baseProfile)
        {
            switch (baseProfile)
            {
                case "T5.5": return "T8.5";
                case "T7.5": return "T10.5";
                case "T8.5": return "T12.5";
                default: return "T14";
            }
        }

        // --------------------------------------------------------------------
        // SOLID CLEANUP
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

            if (!b.IsSolid)
                info?.Add($"{label}: still not a closed solid after capping.");

            return b;
        }

        // --------------------------------------------------------------------
        // ICON & GUID
        // --------------------------------------------------------------------

        protected override Bitmap Icon
        {
            get
            {
                return new Bitmap(24, 24);
            }
        }

        public override Guid ComponentGuid =>
            new Guid("3E366D3F-9F1E-4A65-B3A0-6C6B9C9E8C42");
    }
}
