using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;

namespace Moria.ReFac
{
    internal static class GhParamUtils
    {
        public static void AddNumberInput(GH_ComponentParamServer ps, string name, string nick, string desc, double def)
        {
            var p = new Param_Number
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.item
            };
            p.SetPersistentData(new GH_Number(def));
            ps.RegisterInputParam(p);
        }

        public static void AddOptionalCurveInput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_Curve
            {
                Name = name,
                NickName = nick,
                Description = desc,
                Access = GH_ParamAccess.item,
                Optional = true
            };
            ps.RegisterInputParam(p);
        }

        public static void AddPointOutput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_Point { Name = name, NickName = nick, Description = desc, Access = GH_ParamAccess.list };
            ps.RegisterOutputParam(p);
        }

        public static void AddCurveOutput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_Curve { Name = name, NickName = nick, Description = desc, Access = GH_ParamAccess.list };
            ps.RegisterOutputParam(p);
        }

        public static void AddBrepOutput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_Brep { Name = name, NickName = nick, Description = desc, Access = GH_ParamAccess.list };
            ps.RegisterOutputParam(p);
        }

        public static void AddTextOutput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_String { Name = name, NickName = nick, Description = desc, Access = GH_ParamAccess.list };
            ps.RegisterOutputParam(p);
        }

        public static void AddGeometryOutput(GH_ComponentParamServer ps, string name, string nick, string desc)
        {
            var p = new Param_Geometry { Name = name, NickName = nick, Description = desc, Access = GH_ParamAccess.list };
            ps.RegisterOutputParam(p);
        }

        public static string SanitizeNick(string s)
        {
            s = (s ?? "").Trim();
            if (s.Length == 0) s = "Udefinert";

            var arr = new char[s.Length];
            for (int i = 0; i < s.Length; i++)
                arr[i] = char.IsLetterOrDigit(s[i]) ? s[i] : '_';

            var nick = new string(arr);
            if (char.IsDigit(nick[0])) nick = "T_" + nick;
            if (nick.Length > 24) nick = nick.Substring(0, 24);
            return nick;
        }

        public static string MakeUniqueNick(string baseNick, HashSet<string> used)
        {
            var nick = baseNick;
            int i = 2;
            while (used.Contains(nick))
            {
                var suffix = "_" + i;
                int cut = Math.Max(1, 24 - suffix.Length);
                nick = (baseNick.Length > cut ? baseNick.Substring(0, cut) : baseNick) + suffix;
                i++;
            }
            used.Add(nick);
            return nick;
        }
    }
}
