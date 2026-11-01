const fs = require('fs');
const content = fs.readFileSync('src/App.tsx', 'utf8');

const targetStart = "                } else {\n                  // Whole protein chain view centered around CoM (0,0,0)";
const targetEnd = "                }\n              })()}\n            </svg>\n          </div>";

const startIndex = content.indexOf(targetStart);
const endIndex = content.indexOf(targetEnd) + targetEnd.length;

if (startIndex === -1 || endIndex === -1 || endIndex < startIndex) {
  console.log('Error finding target block');
  process.exit(1);
}

const replacement = `                } else {
                  // Whole protein chain view centered around CoM (0,0,0) or Compare view
                  const views = viewMode === 'compare' 
                    ? [
                        { label: 'DETERMINISTIC ENGINE', xOffset: -200, variant: 'det' },
                        { label: 'ALPHAFOLD PDB REF', xOffset: 200, variant: 'af' }
                      ]
                    : [
                        { label: '', xOffset: 0, variant: 'det' }
                      ];

                  const focusedIdx = playbackNodeIndex !== null ? playbackNodeIndex : activeIndex;

                  return (
                    <g>
                      {views.map((v, vIdx) => {
                        const localCX = cX + v.xOffset;
                        
                        const projectedChain = centeredChain.map((c, i) => {
                          const perturb = v.variant === 'af' ? Math.sin(i * 0.8) * 45 : 0;
                          const proj = project3DTo2D(c.x + perturb, c.y + (v.variant === 'af' ? Math.cos(i) * 35 : 0), c.z + perturb, localCX, cY, v.variant === 'af' ? 0.85 : 0.9);
                          return { ...c, ...proj };
                        });

                        return (
                          <g key={'view-' + vIdx}>
                            {v.label && (
                              <text x={localCX} y={cY - 220} className="fill-zinc-500 font-mono text-[14px] font-extrabold text-center tracking-widest" textAnchor="middle">
                                {v.label}
                                {v.variant === 'af' && (
                                  <tspan x={localCX} y={cY - 200} className="fill-emerald-500 text-[10px]">LDDT SCORE: 89.4% (GOOD)</tspan>
                                )}
                              </text>
                            )}
                            
                            {/* Axis System */}
                            {(() => {
                              const p0 = project3DTo2D(0, 0, 0, localCX, cY, 0.9);
                              const pX = project3DTo2D(450, 0, 0, localCX, cY, 0.9);
                              const pY = project3DTo2D(0, 450, 0, localCX, cY, 0.9);
                              const pZ = project3DTo2D(0, 0, 450, localCX, cY, 0.9);
                              return (
                                <g opacity="0.9" className="pointer-events-none select-none">
                                  <line x1={p0.x} y1={p0.y} x2={pX.x} y2={pX.y} stroke="#ef4444" strokeWidth="2.5" strokeDasharray="4,2" />
                                  <circle cx={pX.x} cy={pX.y} r="3" fill="#ef4444" />
                                  <text x={pX.x + 6} y={pX.y + 3} className="fill-red-400 font-mono text-[10px] font-extrabold" style={{ textShadow: '0 0 4px rgba(0,0,0,0.8)' }}>X</text>
                                  
                                  <line x1={p0.x} y1={p0.y} x2={pY.x} y2={pY.y} stroke="#10b981" strokeWidth="2.5" strokeDasharray="4,2" />
                                  <circle cx={pY.x} cy={pY.y} r="3" fill="#10b981" />
                                  <text x={pY.x + 6} y={pY.y + 3} className="fill-emerald-400 font-mono text-[10px] font-extrabold" style={{ textShadow: '0 0 4px rgba(0,0,0,0.8)' }}>Y</text>
                                  
                                  <line x1={p0.x} y1={p0.y} x2={pZ.x} y2={pZ.y} stroke="#3b82f6" strokeWidth="2.5" strokeDasharray="4,2" />
                                  <circle cx={pZ.x} cy={pZ.y} r="3" fill="#3b82f6" />
                                  <text x={pZ.x + 6} y={pZ.y + 3} className="fill-blue-400 font-mono text-[10px] font-extrabold" style={{ textShadow: '0 0 4px rgba(0,0,0,0.8)' }}>Z</text>

                                  <circle cx={p0.x} cy={p0.y} r={7} fill="none" stroke="#fbbf24" strokeWidth="1.5" />
                                  <text x={p0.x + 8} y={p0.y - 8} className="fill-yellow-500 font-mono text-[7px]" opacity="0.9">CENTER OF MASS CoM (0,0,0)</text>
                                </g>
                              );
                            })()}

                            {/* Connection links */}
                            {projectedChain.map((p, idx) => {
                              if (idx === 0) return null;
                              const prev = projectedChain[idx - 1];
                              const isLinkFocused = focusedIdx === idx || focusedIdx === idx - 1;
                              const opacityVal = xrayFocus ? (isLinkFocused ? 1.0 : 0.15) : 1.0;
                              const color = v.variant === 'af' ? 'rgba(16, 185, 129, 0.5)' : 'rgba(99, 102, 241, 0.65)';
                              const beadColor = v.variant === 'af' ? '#10b981' : '#f43f5e';

                              return (
                                <g key={'lnk-'+idx} opacity={opacityVal} className="transition-opacity duration-300">
                                  <line
                                    x1={prev.x}
                                    y1={prev.y}
                                    x2={p.x}
                                    y2={p.y}
                                    stroke={color}
                                    strokeWidth="3.5"
                                    strokeLinecap="round"
                                  />
                                  <circle cx={(prev.x + p.x)/2} cy={(prev.y + p.y)/2} r="3" fill={beadColor} />
                                </g>
                              );
                            })}

                            {/* Golden Disulfide covalent paths */}
                            {(() => {
                              const pairs: number[][] = [];
                              projectedChain.forEach((p, idx) => {
                                if (p.node.aminoAcid.id === 'CYS' || p.node.aminoAcid.id === 'SEC') {
                                  projectedChain.forEach((p2, idx2) => {
                                    if (idx2 > idx && (p2.node.aminoAcid.id === 'CYS' || p2.node.aminoAcid.id === 'SEC')) {
                                      pairs.push([idx, idx2]);
                                    }
                                  });
                                }
                              });

                              return pairs.map(([i1, i2], pIdx) => {
                                const p1 = projectedChain[i1];
                                const p2 = projectedChain[i2];
                                const isBridgeFocused = focusedIdx === i1 || focusedIdx === i2;
                                const baseOpacity = proteinVariant === 'patient' ? 0.08 : 1.0;
                                const opacityVal = xrayFocus ? (isBridgeFocused ? 1.0 : 0.12 * baseOpacity) : baseOpacity;

                                return (
                                  <path
                                    key={'ss-bridge-'+pIdx}
                                    d={\`M \${p1.x} \${p1.y} Q \${(p1.x+p2.x)/2 + 25} \${(p1.y+p2.y)/2 - 25} \${p2.x} \${p2.y}\`}
                                    stroke="#eab308"
                                    strokeWidth="3"
                                    strokeDasharray="4,2"
                                    fill="none"
                                    className="animate-pulse"
                                    opacity={opacityVal}
                                  />
                                );
                              });
                            })()}

                            {/* Ball nodes for each residue */}
                            {projectedChain.map((p, idx) => {
                              const node = p.node;
                              const radius = 17 * p.scale;
                              const isActive = idx === activeIndex;
                              const isPlaying = playbackNodeIndex === idx;
                              const isNodeFocused = idx === focusedIdx;
                              const opacityVal = xrayFocus ? (isNodeFocused ? 1.0 : 0.22) : 1.0;
                              const atomColor = v.variant === 'af' ? '#064e3b' : node.aminoAcid.color;
                              const atomStroke = v.variant === 'af' ? '#10b981' : (isActive ? '#ffffff' : '#000000');
                              const textStyle = v.variant === 'af' ? 'fill-emerald-400' : 'fill-white';

                              return (
                                <g 
                                  key={'chain-nd-'+idx} 
                                  className="cursor-pointer transition-opacity duration-300" 
                                  onClick={() => selectAndSonifyNode(idx)}
                                  opacity={opacityVal}
                                >
                                  <circle cx={p.x} cy={p.y} r={radius + 8} fill={atomColor} opacity={isActive ? 0.35 : isPlaying ? 0.25 : 0} />
                                  
                                  {isPlaying && (
                                    <g>
                                      <circle 
                                        cx={p.x} 
                                        cy={p.y} 
                                        r={radius + 14} 
                                        fill="none" 
                                        stroke="#ef4444" 
                                        strokeWidth="2.5" 
                                        className="animate-pulse" 
                                        opacity="0.9" 
                                      />
                                      <circle 
                                        cx={p.x} 
                                        cy={p.y} 
                                        r={radius + 7} 
                                        fill="none" 
                                        stroke="#f43f5e" 
                                        strokeWidth="1.8" 
                                        opacity="0.8" 
                                      />
                                      <circle 
                                        cx={p.x + radius} 
                                        cy={p.y - radius} 
                                        r="5" 
                                        fill="#ef4444" 
                                        stroke="#ffffff" 
                                        strokeWidth="1.2" 
                                      />
                                    </g>
                                  )}

                                  <circle 
                                    cx={p.x} 
                                    cy={p.y} 
                                    r={radius} 
                                    fill={atomColor} 
                                    stroke={atomStroke} 
                                    strokeWidth={isActive ? '2.5' : '1.5'} 
                                  />
                                  <text
                                    x={p.x}
                                    y={p.y + 3.5}
                                    className={\`\${v.variant === 'af' ? 'fill-emerald-400' : 'fill-black'} font-extrabold text-[10px] font-mono text-center select-none pointer-events-none\`}
                                    textAnchor="middle"
                                  >
                                    {node.aminoAcid.code1}
                                  </text>
                                  
                                  <text
                                    x={p.x}
                                    y={p.y - radius - 6}
                                    className={\`\${textStyle} font-mono text-[9px] text-center select-none pointer-events-none font-bold\`}
                                    textAnchor="middle"
                                    opacity={isActive ? 1 : 0.6}
                                    style={{ textShadow: '1px 1px 2px #000' }}
                                  >
                                    {idx+1}.{node.aminoAcid.code3}
                                  </text>
                                </g>
                              );
                            })}
                          </g>
                        );
                      })}
                    </g>
                  );
                }
              })()}
            </svg>
          </div>`;

const newContent = content.substring(0, startIndex) + replacement + content.substring(endIndex);
fs.writeFileSync('src/App.tsx', newContent);
console.log('Success');
