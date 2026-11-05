import { platforma } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import TitrationCurvesPage from "./pages/TitrationCurvesPage.vue";
import KDDistributionPage from "./pages/KDDistributionPage.vue";
import AffinityVsFitPage from "./pages/AffinityVsFitPage.vue";
import NtVariantsPage from "./pages/NtVariantsPage.vue";
import TablePage from "./pages/TablePage.vue";
import VsParentPage from "./pages/VsParentPage.vue";
import SynonymousNoisePage from "./pages/SynonymousNoisePage.vue";
import ReplicateNoisePage from "./pages/ReplicateNoisePage.vue";
import TitrationEdgesPage from "./pages/TitrationEdgesPage.vue";
import ReadDepthPage from "./pages/ReadDepthPage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  app.model.data.customBlockLabel ??= "";

  return {
    progress: () => app.model.outputs.isRunning,
    routes: {
      "/": () => TablePage,
      "/nt-variants": () => NtVariantsPage,
      "/vs-parent": () => VsParentPage,
      "/synonymous-noise": () => SynonymousNoisePage,
      "/replicate-noise": () => ReplicateNoisePage,
      "/titration-edges": () => TitrationEdgesPage,
      "/read-depth": () => ReadDepthPage,
      "/titration-curves": () => TitrationCurvesPage,
      "/kd-distribution": () => KDDistributionPage,
      "/affinity-vs-fit": () => AffinityVsFitPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
