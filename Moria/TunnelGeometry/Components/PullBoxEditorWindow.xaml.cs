using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;

namespace Moria.TunnelGeometry.Components
{
    public partial class PullBoxEditorWindow : Window
    {
        private ObservableCollection<PullBoxInstance> _items;

        public event Action<List<PullBoxInstance>> ApplyRequested;

        public PullBoxEditorWindow(List<PullBoxInstance> instances)
        {
            InitializeComponent();

            _items = new ObservableCollection<PullBoxInstance>(
                (instances ?? new List<PullBoxInstance>()).Select(x => x.Clone()).ToList());

            GridBoxes.ItemsSource = _items;
        }

        private PullBoxInstance Selected => GridBoxes.SelectedItem as PullBoxInstance;

        private List<PullBoxInstance> GetInstancesSorted() =>
            _items.OrderBy(x => x.EffectiveStation).ToList();

        private void CommitEditsHard()
        {
            Keyboard.ClearFocus();
            GridBoxes.CommitEdit(DataGridEditingUnit.Cell, true);
            GridBoxes.CommitEdit(DataGridEditingUnit.Row, true);
            Dispatcher.Invoke(() => { }, DispatcherPriority.Background);
        }

        private void OnApply(object sender, RoutedEventArgs e)
        {
            CommitEditsHard();
            ApplyRequested?.Invoke(GetInstancesSorted());
        }

        private void OnAdd(object sender, RoutedEventArgs e)
        {
            double s = _items.Count > 0 ? _items.Max(x => x.EffectiveStation) + 25.0 : 0.0;
            _items.Add(new PullBoxInstance { Station = s });
        }

        private void OnRemove(object sender, RoutedEventArgs e)
        {
            var sel = Selected;
            if (sel == null) return;
            _items.Remove(sel);
        }

        private void OnClose(object sender, RoutedEventArgs e) => Close();

        private void OnCurrentCellChanged(object sender, EventArgs e)
        {
            try
            {
                GridBoxes.CommitEdit(DataGridEditingUnit.Cell, true);
                GridBoxes.CommitEdit(DataGridEditingUnit.Row, true);
            }
            catch { }
        }

        private void OnCellEditEnding(object sender, DataGridCellEditEndingEventArgs e)
        {
            if (e.EditAction != DataGridEditAction.Commit)
                return;

            if (e.Row?.Item is not PullBoxInstance inst)
                return;

            string header = e.Column?.Header?.ToString() ?? "";

            if (header == "Station" || header == "dStation")
                inst.FollowInterval = false;

            if (header == "L (Z)" || header == "W (X)" || header == "H (Y)" ||
                header == "LongOpenW (Z)" || header == "LongOpenH (0→Y)" ||
                header == "ShortOpenW (X)" || header == "ShortOpenH (0→Y)" ||
                header == "Thk" ||
                header == "OffsetZ" || header == "OffsetXY")
            {
                inst.UseDefaults = false;
            }

            Dispatcher.BeginInvoke(new Action(() =>
            {
                GridBoxes.Items.Refresh();
            }), DispatcherPriority.Background);
        }
    }
}