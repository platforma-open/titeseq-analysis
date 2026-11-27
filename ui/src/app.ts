import { platforma } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { plRefsEqual } from "@platforma-sdk/model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import { computed, watch } from "vue";
import AffinityVsFitPage from "./pages/AffinityVsFitPage.vue";
import RunQualityPage from "./pages/RunQualityPage.vue";
import TablePage from "./pages/TablePage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  app.model.data.customBlockLabel ??= "";

  syncSubtitleLabel(app.model);

  const routes = {
    "/": () => TablePage,
    "/affinity-vs-fit": () => AffinityVsFitPage,
    "/run-quality": () => RunQualityPage,
  };

  // A project saved on a page the block does not have opens on the main page.
  const savedPath = new URL(app.snapshot.navigationState.href, "http://dummy").pathname;
  if (!(savedPath in routes)) void app.navigateTo("/");

  return {
    progress: () => app.model.outputs.isRunning,
    routes,
  };
});

export const useApp = sdkPlugin.useApp;

type AppModel = ReturnType<typeof useApp>["model"];

// The subtitle names the picked abundance as the result pool labels it; the subset's label when
// one is picked, since it already carries the dataset as a prefix. `subtitleLabel` never reaches
// args, so following the pool here cannot stale the block.
function syncSubtitleLabel(model: AppModel) {
  // "" when no abundance is picked; undefined while the picked abundance or subset is not among
  // the options yet.
  const subtitleLabel = computed<string | undefined>(() => {
    const { abundanceRef, filterRef } = model.data;
    if (!abundanceRef) return "";
    const option = model.outputs.abundanceOptions?.find((o) =>
      plRefsEqual(o.primary.ref, abundanceRef, true),
    );
    if (!option) return undefined;
    if (!filterRef) return option.primary.label;
    return option.filters?.find((f) => plRefsEqual(f.ref, filterRef, true))?.label;
  });

  watch(
    subtitleLabel,
    (label) => {
      if (label !== undefined && label !== model.data.subtitleLabel)
        model.data.subtitleLabel = label;
    },
    { immediate: true },
  );
}
