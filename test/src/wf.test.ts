import type {
  BlockOutputs,
  model,
} from '@platforma-open/platforma-open.titeseq-analysis.model';
import type { InferBlockState } from '@platforma-sdk/model';
import { wrapOutputs } from '@platforma-sdk/model';
import { awaitStableState, blockTest } from '@platforma-sdk/test';
import { blockSpec as myBlockSpec } from 'this-block';

blockTest('empty inputs', { timeout: 20000 }, async ({ rawPrj: project, expect }) => {
  const blockId = await project.addBlock('Titeseq Analysis', myBlockSpec);

  const stableState = (await awaitStableState(
    project.getBlockState(blockId),
    15000,
  )) as InferBlockState<typeof model>;

  expect(stableState.outputs).toMatchObject({
    abundanceOptions: { ok: true, value: [] },
    concentrationOptions: { ok: true, value: [] },
    binOptions: { ok: true, value: [] },
    antigenOptions: { ok: true, value: [] },
    validationWarnings: { ok: true, value: [] },
    binMode: { ok: true, value: false },
  });

  const outputs = wrapOutputs<BlockOutputs>(stableState.outputs);
  expect(outputs.abundanceOptions).toEqual([]);
  expect(outputs.concentrationOptions).toEqual([]);
  expect(outputs.binOptions).toEqual([]);
  expect(outputs.antigenOptions).toEqual([]);
  expect(outputs.validationWarnings).toEqual([]);
  expect(outputs.binMode).toEqual(false);

  // Optional workflow outputs must not appear before the block runs (argsValid is false).
  expect(stableState.outputs).not.toHaveProperty('summaryPfHandle');
  expect(stableState.outputs).not.toHaveProperty('titrationCurvesPf');
  expect(stableState.outputs).not.toHaveProperty('logHandle');

  // isRunning should be falsy (undefined or false) while the block is idle.
  expect(outputs.isRunning).toBeFalsy();
});
