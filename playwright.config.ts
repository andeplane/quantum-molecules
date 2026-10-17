import { defineConfig } from '@playwright/test';

// The app is opened with `?headless`: headless Chromium runs WebGPU on SwiftShader but cannot
// present a WebGPU canvas, so the tests render offscreen and read the pixels back.
export default defineConfig({
  testDir: 'e2e',
  timeout: 120_000,
  fullyParallel: false,
  workers: 1,
  reporter: process.env.CI ? 'list' : 'line',
  use: {
    baseURL: 'http://localhost:4173/',
    viewport: { width: 1280, height: 800 },
    launchOptions: {
      args: ['--enable-unsafe-webgpu', '--enable-unsafe-swiftshader', '--use-webgpu-adapter=swiftshader'],
    },
  },
  webServer: {
    command: 'npx vite build && npx vite preview --port 4173 --strictPort',
    url: 'http://localhost:4173/',
    reuseExistingServer: !process.env.CI,
    timeout: 120_000,
  },
});
