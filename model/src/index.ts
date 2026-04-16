import type { InferOutputsType } from "@platforma-sdk/model";
import { BlockModelV3, DataModelBuilder } from "@platforma-sdk/model";

export type BlockData = {
  name: string;
};

const dataModel = new DataModelBuilder().from<BlockData>("v1").init(() => ({ name: "" }));

export const platforma = BlockModelV3.create(dataModel)

  .args((data) => ({ name: data.name }))

  .output("tengoMessage", (ctx) => ctx.outputs?.resolve("tengoMessage")?.getDataAsJson())

  .output("pythonMessage", (ctx) => ctx.outputs?.resolve("pythonMessage")?.getDataAsString())

  .sections((_ctx) => [{ type: "link", href: "/", label: "Main" }])

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
