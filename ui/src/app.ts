import { platforma } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import ExplorerPage from "./pages/ExplorerPage.vue";
import NtVariantsPage from "./pages/NtVariantsPage.vue";
import QcPage from "./pages/QcPage.vue";
import ReplicatesPage from "./pages/ReplicatesPage.vue";
import TablePage from "./pages/TablePage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  app.model.data.customBlockLabel ??= "";

  return {
    progress: () => app.model.outputs.isRunning,
    routes: {
      "/": () => TablePage,
      "/nt-variants": () => NtVariantsPage,
      "/explorer": () => ExplorerPage,
      "/replicates": () => ReplicatesPage,
      "/qc": () => QcPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
