using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace Moria.ReFac
{
    public class MoriaInfo : GH_AssemblyInfo
    {
        public override string Name => "Moria";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("5b327509-d34d-43c2-a85e-19df1e548d8f");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";

        //Return a string representing the version.  This returns the same version as the assembly.
        public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
    }
}