const fs = require('fs');
let content = fs.readFileSync('src/App.tsx', 'utf8');

const replacements = [
  ['Pitch (Oś X):', 'Pitch (X Axis):'],
  ['Yaw (Oś Y):', 'Yaw (Y Axis):'],
  ['Myszka rotuje / Scroll przybliża', 'Drag to rotate / Scroll to zoom'],
  ['OGNIWO NR', 'RESIDUE #'],
  ['Masa molowa:', 'Molar Mass:'],
  ['Wzór:', 'Formula:'],
  ['STAN: ZABLOKOWANE (SNAPLOCK)', 'STATE: LOCKED (SNAPLOCK)'],
  ['STAN: GIĘTKIE (LOSOWANIA)', 'STATE: FLEXIBLE (KINETICS)'],
  ['Model białka zablokowany sferycznie względem ciężkości CoM.', 'Protein model spherically locked to Center of Mass.'],
  ['Kamera:', 'Camera:'],
  ['Przywróć domyślne ustawienia kamery', 'Restore default camera settings'],
  ['Skala (Zoom):', 'Scale (Zoom):'],
  ['Oddal (-)', 'Zoom Out (-)'],
  ['Zwiń deterministycznie', 'Fold Deterministically'],
  ['Zwijanie fizyczne...', 'Folding kinetics simulation...'],
];

for (const [pl, en] of replacements) {
  content = content.replace(pl, en);
}

fs.writeFileSync('src/App.tsx', content);
console.log('Done');
