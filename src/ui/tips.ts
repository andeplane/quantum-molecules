/** Things a newcomer does, in the order the tips suggest them. */
export type TipEvent = 'drag' | 'relax' | 'orbital' | 'preset';

export interface Tip {
  /** The action that completes this tip. */
  until: TipEvent;
  /** HTML; written here, never from user input. */
  html: string;
}

export const TIPS: Tip[] = [
  { until: 'drag', html: '<b>Try it:</b> drag one of the atoms away from the others. Every frame is a fresh quantum calculation.' },
  {
    until: 'relax',
    html: 'The energy on the right went up: the molecule is strained. Press <b>Relax geometry</b> and watch the atoms settle back.',
  },
  { until: 'orbital', html: 'Now look inside. Click a line under <b>Energy levels</b> on the right to see where that pair of electrons lives.' },
  {
    until: 'preset',
    html: 'That shape is an orbital: blue and orange are the two signs of its wave. Next, pick another molecule from the list on the left.',
  },
];

/** The first tip whose action has not happened yet, or null when all are done. */
export function currentTip(done: ReadonlySet<TipEvent>): Tip | null {
  return TIPS.find((t) => !done.has(t.until)) ?? null;
}
