import { fileURLToPath } from 'node:url';
import { defineConfig } from 'vitest/config';

const qchemPkg = fileURLToPath(new URL('./crates/qchem-wasm/pkg', import.meta.url));

export default defineConfig({
  base: './',
  resolve: { alias: { '@qchem': qchemPkg } },
  worker: { format: 'es' },
  test: {
    include: ['tests/**/*.test.ts'],
    testTimeout: 120000,
  },
});
