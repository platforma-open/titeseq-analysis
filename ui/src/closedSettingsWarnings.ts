import { ref } from "vue";

/**
 * Settings warnings the user closed. Module scope, so a warning stays closed when the drawer closes
 * and opens again, until the block's page reloads.
 */
export const closedSettingsWarnings = ref<ReadonlySet<string>>(new Set());
