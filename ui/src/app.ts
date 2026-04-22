import { model } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import TitrationCurvesPage from "./pages/TitrationCurvesPage.vue";
import KDDistributionPage from "./pages/KDDistributionPage.vue";
import AffinityVsFitPage from "./pages/AffinityVsFitPage.vue";
import TablePage from "./pages/TablePage.vue";

export const sdkPlugin = defineAppV3(model, (app) => {
  app.model.data.customBlockLabel ??= "";

  return {
    progress: () => app.model.outputs.isRunning,
    routes: {
      "/": () => TablePage,
      "/titration-curves": () => TitrationCurvesPage,
      "/kd-distribution": () => KDDistributionPage,
      "/affinity-vs-fit": () => AffinityVsFitPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
