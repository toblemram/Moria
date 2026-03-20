using System;
using System.Collections.Generic;
using System.Linq;

namespace Moria.ReFac
{
    internal static class SosiUtils
    {
        public static string NormalizeType(string objType)
        {
            var t = (objType ?? "").Trim();
            return t.Length == 0 ? "Udefinert" : t;
        }

        public static List<string> CollectAllObjTypes(SosiData d)
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

            if (set.Count == 0) anyUndef = true;

            var list = set.OrderBy(s => s, StringComparer.OrdinalIgnoreCase).ToList();
            if (anyUndef) list.Insert(0, "Udefinert");

            return list;
        }
    }
}
