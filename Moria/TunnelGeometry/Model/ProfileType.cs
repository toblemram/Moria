using System.Collections.Generic;

namespace Moria.TunnelGeometry.Components
{
    /// <summary>
    /// T-profile parameters as given in the handbook table.
    /// All lengths in metres.
    /// </summary>
    public static class ProfileType
    {
        public readonly struct ProfileParameters
        {
            public double Yv { get; }
            public double Rv { get; }
            public double X { get; }
            public double Rh { get; }

            public ProfileParameters(double yv, double rv, double x, double rh)
            {
                Yv = yv;
                Rv = rv;
                X = x;
                Rh = rh;
            }
        }

        /// <summary>
        /// Profiles: key = "T5.5", "T7.5", ..., "T14".
        /// Value = (Yv, Rv, X, Rh).
        /// </summary>
        public static readonly IReadOnlyDictionary<string, ProfileParameters> Profiles =
            new Dictionary<string, ProfileParameters>
            {
                // LOW-ROOF PROFILES (special construction)
                { "T5.5",  new ProfileParameters(1.770, 4.790, 3.402, 2.587) },
                { "T7.5",  new ProfileParameters(1.570, 4.790, 1.550, 3.594) },
                { "T8.5",  new ProfileParameters(1.770, 4.790, 0.402, 4.500) },

                // STANDARD PROFILES (same construction as old component)
                { "T9.5",  new ProfileParameters(1.570, 4.790, 0.450, 5.212) },
                { "T10.5", new ProfileParameters(1.570, 4.790, 1.450, 5.950) },
                { "T11.5", new ProfileParameters(1.770, 4.790, 2.598, 7.199) },
                { "T12.5", new ProfileParameters(1.570, 4.790, 3.450, 7.458) },
                { "T13",   new ProfileParameters(1.570, 4.790, 3.950, 7.825) },
                { "T13.5", new ProfileParameters(1.570, 4.790, 4.450, 8.053) },
                { "T14",   new ProfileParameters(1.570, 4.790, 4.950, 8.575) },
            };

        /// <summary>
        /// Tabulated Yh for low-roof profiles (centre height of Rh).
        /// These do NOT follow the same formula as the big profiles.
        /// </summary>
        private static readonly Dictionary<string, double> LowRoofYh =
            new Dictionary<string, double>
            {
                { "T5.5", 3.171 },
                { "T7.5", 2.481 },
                { "T8.5", 1.981 },
            };

        public static bool IsLowRoof(string type) =>
            LowRoofYh.ContainsKey(type);

        public static bool TryGetLowRoofYh(string type, out double yh) =>
            LowRoofYh.TryGetValue(type, out yh);
    }
}
