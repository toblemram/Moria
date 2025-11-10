using System;
using System.Text;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Tunnel.GH
{
    /// <summary>
    /// Komponent som forsøker robust å lukke en tunnel-mesh
    /// (HealNakedEdges, FillHoles, ShrinkWrap), og i tillegg
    /// forsøker å lage solider:
    ///  - MeshSolid: Brep fra det lukkede meshet
    ///  - BrepSolid: Inngangsbrep gjort så solid som mulig
    /// </summary>
    public class GH_TunnelMeshCloser : GH_Component
    {
        public GH_TunnelMeshCloser()
          : base(
                "Tunnel.MeshCloser",
                "MeshCloser",
                "Forsøker robust å lukke en tunnelmesh (HealNakedEdges, FillHoles, ShrinkWrap), og lager solider for både mesh og brep.",
                "Tunnel",
                "Analyse")
        { }

        public override Guid ComponentGuid => new Guid("F2BBD9E8-3D0A-4B4D-9B0B-7C1C0E8B7C21");
        protected override System.Drawing.Bitmap Icon => null;

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddMeshParameter("Mesh", "M", "Mesh som skal forsøkes lukket.", GH_ParamAccess.item);
            p.AddBrepParameter("DesignBrep", "B", "Brep (prosjektert tunnel el.l.) som forsøkes gjort solid.", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddMeshParameter("Closed Mesh", "CM", "Resultatmesh etter lukke-forsøk.", GH_ParamAccess.item);
            p.AddBrepParameter("Mesh Solid", "MS", "Solid Brep laget fra det (tilnærmet) lukkede meshet.", GH_ParamAccess.item);
            p.AddBrepParameter("Brep Solid", "BS", "DesignBrep forsøkt gjort solid (CapPlanarHoles).", GH_ParamAccess.item);
            p.AddBooleanParameter("Is Closed", "C", "True hvis sluttmesh er lukket.", GH_ParamAccess.item);
            p.AddTextParameter("Report", "R", "Logg over hva som ble forsøkt og resultat.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            Mesh inputMesh = null;
            Brep inputBrep = null;

            if (!da.GetData(0, ref inputMesh) || inputMesh == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Ingen mesh gitt inn.");
                return;
            }
            if (!da.GetData(1, ref inputBrep) || inputBrep == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Ingen brep gitt inn.");
                return;
            }

            var log = new StringBuilder();

            // Toleranse
            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3;
            log.AppendLine($"Tol = {tol:0.###}");

            // Kopi så vi ikke ødelegger originalen
            Mesh mesh = inputMesh.DuplicateMesh();

            log.AppendLine("=== Tunnel mesh lukker ===");
            log.AppendLine($"Start: Vertices: {mesh.Vertices.Count}, Faces: {mesh.Faces.Count}");

            // 1) Sjekk gyldighet
            string validLog;
            bool isValid = mesh.IsValidWithLog(out validLog);
            if (!isValid)
            {
                log.AppendLine("Mesh er IKKE gyldig. Forsøker enkel reparasjon...");
                log.AppendLine(validLog);

                // Enkle reparasjoner
                mesh.Vertices.CullUnused();
                mesh.Weld(Math.PI); // sveiser sammen verts som deler posisjon
                mesh.UnifyNormals();
                mesh.Normals.ComputeNormals();
                mesh.FaceNormals.ComputeFaceNormals();
                mesh.Compact();

                string validLog2;
                bool isValid2 = mesh.IsValidWithLog(out validLog2);
                if (!isValid2)
                {
                    log.AppendLine("Klarte ikke å gjøre mesh gyldig automatisk.");
                    log.AppendLine(validLog2);

                    // Output med minimal behandling
                    da.SetData(0, mesh);
                    da.SetData(1, null);  // MeshSolid
                    da.SetData(2, TryMakeSolidFromBrep(inputBrep, tol, log, "BrepSolid"));
                    da.SetData(3, false);
                    da.SetData(4, log.ToString());
                    return;
                }
                else
                {
                    log.AppendLine("Lyktes i å fikse noen mesh-feil automatisk.");
                }
            }
            else
            {
                log.AppendLine("Mesh er gyldig.");
            }

            // Hjelpefunksjon for å telle naked edges
            Func<Mesh, int> countNakedEdges = (mm) =>
            {
                var naked = mm.GetNakedEdges();
                if (naked == null) return 0;
                int c = 0;
                foreach (var pl in naked)
                    c += Math.Max(0, pl.Count - 1);
                return c;
            };

            int nakedStart = countNakedEdges(mesh);
            log.AppendLine($"Naked edges ved start (ca): {nakedStart}");

            if (mesh.IsClosed && nakedStart == 0)
            {
                log.AppendLine("Mesh var allerede lukket.");
            }
            else
            {
                // 2) HealNakedEdges – sveise nesten-sammen kantene
                BoundingBox bb = mesh.GetBoundingBox(true);
                double diag = bb.Diagonal.Length;
                double healDist = diag * 0.001; // 0.1% av modellstørrelse – juster ved behov

                try
                {
                    bool healed = mesh.HealNakedEdges(healDist);
                    log.AppendLine($"HealNakedEdges({healDist:0.###}) => {(healed ? "true" : "false")}");
                }
                catch (Exception e)
                {
                    log.AppendLine("HealNakedEdges feilet: " + e.Message);
                }

                int nakedAfterHeal = countNakedEdges(mesh);
                log.AppendLine($"Naked edges etter HealNakedEdges: {nakedAfterHeal}");

                // 3) FillHoles – fyll hull definert av naked edges
                try
                {
                    bool filled = mesh.FillHoles();
                    log.AppendLine("FillHoles() forsøkt: " + (filled ? "fylte ett eller flere hull." : "ingen hull ble fylt."));
                }
                catch (Exception e)
                {
                    log.AppendLine("FillHoles feilet: " + e.Message);
                }

                int nakedAfterFill = countNakedEdges(mesh);
                log.AppendLine($"Naked edges etter FillHoles: {nakedAfterFill}");

                // 4) ShrinkWrap (Rhino 8+): lag en ny, vanntett mesh rundt input
                Mesh shrinkWrapped = null;
                bool shrinkSuccess = false;

                try
                {
                    // NB: krever Rhino 8 og RhinoCommon med ShrinkWrapParameters
                    var p = new ShrinkWrapParameters
                    {
                        TargetEdgeLength = diag * 0.01,   // ca 1% av diag – juster for mer/mindre detalj
                        FillHolesInInputObjects = true,
                        SmoothingIterations = 3,
                        Offset = 0.0,
                        PolygonOptimization = 0          // maks detalj
                    };

                    shrinkWrapped = mesh.ShrinkWrap(p);
                    if (shrinkWrapped != null)
                    {
                        shrinkSuccess = true;
                        log.AppendLine("ShrinkWrap laget en ny mesh.");
                    }
                    else
                    {
                        log.AppendLine("ShrinkWrap returnerte null (ingen mesh).");
                    }
                }
                catch (Exception e)
                {
                    log.AppendLine("ShrinkWrap feilet (bruker du Rhino 8?): " + e.Message);
                }

                if (shrinkSuccess)
                {
                    mesh = shrinkWrapped;
                    int nakedShrink = countNakedEdges(mesh);
                    log.AppendLine($"Naked edges etter ShrinkWrap: {nakedShrink}");
                }
            }

            // 5) Endelig sjekk – er den lukket?
            int nakedFinal = countNakedEdges(mesh);
            bool closed = mesh.IsClosed && nakedFinal == 0;

            log.AppendLine("=== Sluttstatus for mesh ===");
            log.AppendLine($"IsClosed (Mesh.IsClosed && naked==0): {closed}");
            log.AppendLine($"Naked edges (final): {nakedFinal}");

            string finalValidLog;
            mesh.IsValidWithLog(out finalValidLog);
            log.AppendLine("Mesh-validitetslogg (til slutt):");
            log.AppendLine(finalValidLog);

            if (!closed)
            {
                log.AppendLine();
                log.AppendLine("Klarte ikke å få en helt lukket mesh.");
                log.AppendLine("Vanlige årsaker:");
                log.AppendLine("- Veldig store eller kompliserte hull som er vanskelige å triangulere.");
                log.AppendLine("- Ikke-manifold kanter (kanter delt av mer enn to flater).");
                log.AppendLine("- Flere frakoblede mesh-deler / overlappende geometri.");
                log.AppendLine();
                log.AppendLine("Tips:");
                log.AppendLine("- Bruk _ShowEdges i Rhino (Naked edges) for å se hvor det fortsatt er åpninger.");
                log.AppendLine("- Vurder å rydde bort “løse” biter av skannet mesh før du kjører komponenten.");
                log.AppendLine("- Test evt. annen oppløsning på ShrinkWrap (TargetEdgeLength).");
            }

            // 6) Forsøk å lage solider

            // 6a) Mesh -> Brep (MeshSolid)
            Brep meshSolid = null;
            try
            {
                meshSolid = Brep.CreateFromMesh(mesh, true);
                if (meshSolid == null || !meshSolid.IsValid)
                {
                    log.AppendLine("Mesh→Brep: CreateFromMesh ga null/ugyldig.");
                    meshSolid = null;
                }
                else
                {
                    if (!meshSolid.IsSolid)
                    {
                        var capped = meshSolid.CapPlanarHoles(tol);
                        if (capped != null && capped.IsValid)
                            meshSolid = capped;

                        if (!meshSolid.IsSolid)
                            log.AppendLine("MeshSolid er fortsatt ikke solid etter CapPlanarHoles.");
                        else
                            log.AppendLine("MeshSolid ble solid etter CapPlanarHoles.");
                    }
                    else
                    {
                        log.AppendLine("MeshSolid var solid direkte fra CreateFromMesh.");
                    }
                }
            }
            catch (Exception e)
            {
                log.AppendLine("Mesh→Brep feilet: " + e.Message);
                meshSolid = null;
            }

            // 6b) Brep -> BrepSolid
            Brep brepSolid = TryMakeSolidFromBrep(inputBrep, tol, log, "BrepSolid");

            // 7) Output
            da.SetData(0, mesh);
            da.SetData(1, meshSolid);
            da.SetData(2, brepSolid);
            da.SetData(3, closed);
            da.SetData(4, log.ToString());
        }

        private Brep TryMakeSolidFromBrep(Brep input, double tol, StringBuilder log, string label)
        {
            if (input == null || !input.IsValid)
            {
                log.AppendLine($"{label}: Inngangsbrep er null/ugyldig.");
                return null;
            }

            var b = input.DuplicateBrep();

            if (b.IsSolid)
            {
                log.AppendLine($"{label}: Inngangsbrep er allerede solid.");
                return b;
            }

            try
            {
                var capped = b.CapPlanarHoles(tol);
                if (capped != null && capped.IsValid)
                {
                    b = capped;
                    if (b.IsSolid)
                        log.AppendLine($"{label}: Ble solid etter CapPlanarHoles.");
                    else
                        log.AppendLine($"{label}: CapPlanarHoles ga brep, men den er fortsatt ikke solid.");
                }
                else
                {
                    log.AppendLine($"{label}: CapPlanarHoles ga null/ugyldig.");
                }
            }
            catch (Exception e)
            {
                log.AppendLine($"{label}: CapPlanarHoles feilet: " + e.Message);
            }

            return b;
        }
    }
}
