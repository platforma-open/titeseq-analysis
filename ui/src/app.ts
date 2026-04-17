import { model } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { defineApp } from "@platforma-sdk/ui-vue";
import TitrationCurvesPage from "./pages/TitrationCurvesPage.vue";
import KDDistributionPage from "./pages/KDDistributionPage.vue";
import AffinityVsFitPage from "./pages/AffinityVsFitPage.vue";
import TablePage from "./pages/TablePage.vue";

export const sdkPlugin = defineApp(model, (app) => {
  app.model.args.customBlockLabel ??= "";

  return {
    progress: () => app.model.outputs.isRunning,
    routes: {
      "/": () => TitrationCurvesPage,
      "/kd-distribution": () => KDDistributionPage,
      "/affinity-vs-fit": () => AffinityVsFitPage,
      "/table": () => TablePage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
