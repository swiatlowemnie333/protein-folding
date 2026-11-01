import React, { useState, useRef, useEffect } from 'react';
import { 
  Bot, 
  Send, 
  Sparkles, 
  Download, 
  Copy, 
  Check, 
  Trash2, 
  Cpu, 
  Layers, 
  Atom, 
  AlertCircle,
  Minimize2,
  Maximize2,
  X as CloseIcon,
  ChevronUp
} from 'lucide-react';
import { generateAminoAcidWarehouseJson } from '../snaplockEngine';

interface Message {
  id: string;
  role: 'user' | 'model';
  text: string;
  timestamp: string;
  modelUsed?: string;
}

interface SnaplockChatbotProps {
  onClose?: () => void;
  isEmbedded?: boolean;
}

export const SnaplockChatbot: React.FC<SnaplockChatbotProps> = ({ onClose, isEmbedded = false }) => {
  const [messages, setMessages] = useState<Message[]>([
    {
      id: 'welcome',
      role: 'model',
      text: `Dzień dobry. Jestem **Głównym Inżynierem Biofizyki Obliczeniowej i Architektem Silnika SNAPLOCK**. 

Wszystkie obliczenia opieram na **dyskretnej algebrze Liego E8 (248 wektorów bazowych i 2232 mikrosiatek)** oraz wgranej matrycy **22 parametrów fizykochemicznych** dla 21 aminokwasów (w tym Selenocysteiny SEC/U).

**Główne reguły geometryczne SNAPLOCK:**
- Port wejściowy N (wiązanie peptydowe): **Łatka #1 (+Y)** \`[0, 1, 0]\`
- Port wyjściowy C (karbonyl): **Łatka #248 (-Y)** \`[0, -1, 0]\`
- Orientacja boczna: **#50 (+Z)**, **#100 (+X)**, **#150 (-Z)**, **#200 (-X)**
- Ewaluacja złączy: **5 ścieżek (Tracks 1-5)** z najwyższym priorytetem dla klamer kowalencyjnych S-S i Se-Se.

Możesz kliknąć poniższy przycisk, aby wygenerować plik \`amino_acid_warehouse.json\` dla pierwszych trzech aminokwasów (GLY, ALA, VAL) lub zapytać o dowolny aspekt mechaniki zwijania.`,
      timestamp: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' }),
      modelUsed: 'gemini-3.8-flash'
    }
  ]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [selectedModel, setSelectedModel] = useState<'gemini-3.8-flash' | 'gemini-4-argon' | 'gemini-3.1-flash-lite'>('gemini-3.8-flash');
  const [copiedId, setCopiedId] = useState<string | null>(null);
  const [isMinimized, setIsMinimized] = useState(false);
  const messagesEndRef = useRef<HTMLDivElement>(null);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    if (!isMinimized) {
      scrollToBottom();
    }
  }, [messages, isLoading, isMinimized]);

  const quickPrompts = [
    "Bazując na załadowanej tabeli 22 parametrów, wygeneruj plik amino_acid_warehouse.json dla pierwszych trzech aminokwasów, mapując ich właściwości na siatkę 248 łatek (E8).",
    "Przeprowadź ewaluację 5 ścieżek (Tracks 1-5) dla klamry CYS-CYS (mostek S-S) vs SEC-SEC (mostek Se-Se).",
    "Wyjaśnij geometrię portów #1 (+Y), #248 (-Y) oraz orientację wektorów #50, #100, #150, #200.",
    "Porównaj objętość Van der Waalsa oraz pKr Selenocysteiny (SEC/U) i Cysteiny (CYS/C)."
  ];

  const handleSend = async (textToSend?: string) => {
    const promptText = (textToSend !== undefined ? textToSend : input).trim();
    if (!promptText || isLoading) return;

    const userMsg: Message = {
      id: 'usr_' + Date.now(),
      role: 'user',
      text: promptText,
      timestamp: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })
    };

    const newHistory = [...messages, userMsg];
    setMessages(newHistory);
    setInput('');
    setIsLoading(true);

    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 10000); // 10s timeout

    try {
      const response = await fetch('/api/gemini/chat', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          messages: newHistory.map(m => ({ role: m.role, text: m.text })),
          modelName: selectedModel
        }),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();
      const modelMsg: Message = {
        id: 'bot_' + Date.now(),
        role: 'model',
        text: data.text || 'Otrzymano pustą odpowiedź od silnika SNAPLOCK.',
        timestamp: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' }),
        modelUsed: data.model || selectedModel
      };

      setMessages(prev => [...prev, modelMsg]);
    } catch (err: any) {
      clearTimeout(timeoutId);
      console.warn('Chat API fallback:', err);
      // Fallback response directly from local SNAPLOCK generator
      let fallbackText = `Informacja silnika SNAPLOCK E8: zapytanie obsłużone lokalnie przez deterministyczny moduł biofizyczny.\n\n`;
      if (promptText.toLowerCase().includes('amino_acid_warehouse.json') || promptText.toLowerCase().includes('warehouse')) {
        const warehouse = generateAminoAcidWarehouseJson(['GLY', 'ALA', 'VAL']);
        fallbackText += `Oto wygenerowany plik magazynu \`amino_acid_warehouse.json\` dla GLY, ALA, VAL z mapowaniem na 248 łatek E8:\n\n\`\`\`json\n${JSON.stringify(warehouse, null, 2)}\n\`\`\``;
      } else {
        fallbackText += `Deterministyczny silnik SNAPLOCK potwierdza gotowość geometryczną. Wszystkie złącza peptydowe operują w sztywnym układzie portu #248 (C-koniec) do portu #1 (N-koniec) z dyskretną rotacją wokół osi styku oraz klamrami disiarczkowymi CYS-CYS (5.5 Å Cα-Cα).`;
      }

      setMessages(prev => [
        ...prev,
        {
          id: 'bot_' + Date.now(),
          role: 'model',
          text: fallbackText,
          timestamp: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' }),
          modelUsed: 'snaplock-local'
        }
      ]);
    } finally {
      setIsLoading(false);
    }
  };

  const handleCopy = (text: string, id: string) => {
    navigator.clipboard.writeText(text);
    setCopiedId(id);
    setTimeout(() => setCopiedId(null), 2000);
  };

  const handleDownloadJson = (jsonContent: string, filename = 'amino_acid_warehouse.json') => {
    const blob = new Blob([jsonContent], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const extractJsonFromMessage = (text: string): string | null => {
    const jsonMatch = text.match(/```(?:json)?\s*([\s\S]*?)\s*```/);
    if (jsonMatch && jsonMatch[1]) {
      try {
        JSON.parse(jsonMatch[1]);
        return jsonMatch[1];
      } catch {
        return null;
      }
    }
    return null;
  };

  // IF MINIMIZED: Render compact floating pill at bottom-right that DOES NOT cover the 3D screen!
  if (!isEmbedded && isMinimized) {
    return (
      <div className="fixed bottom-4 right-4 z-50 flex items-center gap-2 bg-slate-950/95 backdrop-blur-md border border-cyan-500/60 rounded-full px-3.5 py-2 shadow-2xl font-mono text-xs text-white">
        <button
          onClick={() => setIsMinimized(false)}
          className="flex items-center gap-2 text-cyan-300 hover:text-white transition cursor-pointer"
          title="Rozwiń czat SNAPLOCK AI"
        >
          <div className="w-6 h-6 rounded-full bg-cyan-950 border border-cyan-500/50 flex items-center justify-center text-cyan-400">
            <Bot className="w-3.5 h-3.5" />
          </div>
          <span className="font-bold">Architekt SNAPLOCK AI</span>
          <ChevronUp className="w-4 h-4 text-cyan-400 animate-bounce" />
        </button>

        {onClose && (
          <button
            onClick={onClose}
            className="p-1 rounded-full text-slate-400 hover:text-rose-400 hover:bg-slate-900 transition-colors cursor-pointer ml-1"
            title="Zamknij czat"
          >
            <CloseIcon className="w-3.5 h-3.5" />
          </button>
        )}
      </div>
    );
  }

  return (
    <div className={`flex flex-col bg-slate-900 border border-slate-800 rounded-2xl shadow-2xl overflow-hidden ${
      isEmbedded ? 'h-full w-full' : 'fixed bottom-4 right-4 w-[520px] max-w-[94vw] h-[620px] max-h-[85vh] z-50'
    }`}>
      {/* Header */}
      <div className="flex items-center justify-between px-4 py-3 bg-slate-950 border-b border-slate-800 shrink-0">
        <div className="flex items-center gap-3">
          <div className="w-9 h-9 rounded-xl bg-gradient-to-tr from-cyan-600 via-indigo-600 to-amber-500 p-0.5 shadow-md flex items-center justify-center">
            <div className="w-full h-full bg-slate-950 rounded-[10px] flex items-center justify-center">
              <Bot className="w-5 h-5 text-cyan-400" />
            </div>
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h3 className="font-sans font-bold text-sm text-slate-100 tracking-wide">
                SNAPLOCK Architect
              </h3>
              <span className="text-[10px] px-1.5 py-0.5 rounded bg-cyan-950/80 text-cyan-300 border border-cyan-800/60 font-mono">
                E8 Lie • 248 Roots
              </span>
            </div>
            <p className="text-[11px] text-slate-400">
              Główny Inżynier Biofizyki Obliczeniowej
            </p>
          </div>
        </div>

        {/* Model selector & window controls */}
        <div className="flex items-center gap-2">
          <div className="relative">
            <select
              value={selectedModel}
              onChange={(e) => setSelectedModel(e.target.value as any)}
              className="text-xs bg-slate-900 text-slate-300 border border-slate-700 rounded-lg px-2 py-1 pr-6 focus:outline-none focus:border-cyan-500 font-mono cursor-pointer"
            >
              <option value="gemini-3.8-flash">gemini-3.8-flash (Standard)</option>
              <option value="gemini-3.1-flash-lite">gemini-3.1-flash-lite (Fast)</option>
              <option value="gemini-4-argon">gemini-4-argon (Advanced QM/MM)</option>
            </select>
          </div>

          {!isEmbedded && (
            <button
              onClick={() => setIsMinimized(true)}
              className="p-1.5 text-slate-400 hover:text-slate-200 hover:bg-slate-800 rounded transition-colors cursor-pointer"
              title="Zwiń do paska"
            >
              <Minimize2 className="w-4 h-4" />
            </button>
          )}

          {onClose && (
            <button
              onClick={onClose}
              className="p-1.5 text-slate-400 hover:text-rose-400 hover:bg-slate-800 rounded transition-colors cursor-pointer"
              title="Zamknij czat"
            >
              <CloseIcon className="w-4 h-4" />
            </button>
          )}
        </div>
      </div>

      {/* Messages scroll area */}
      <div className="flex-1 overflow-y-auto p-4 space-y-4 bg-slate-900/60">
        {messages.map((m) => {
          const detectedJson = m.role === 'model' ? extractJsonFromMessage(m.text) : null;

          return (
            <div
              key={m.id}
              className={`flex flex-col ${m.role === 'user' ? 'items-end' : 'items-start'}`}
            >
              <div className="flex items-center gap-1.5 mb-1 px-1">
                <span className="text-[10px] font-mono font-medium text-slate-400">
                  {m.role === 'user' ? 'Ty' : 'Główny Inżynier SNAPLOCK'}
                </span>
                <span className="text-[9px] text-slate-500">{m.timestamp}</span>
                {m.modelUsed && (
                  <span className="text-[9px] font-mono px-1 py-0.2 rounded bg-slate-800 text-slate-400">
                    {m.modelUsed}
                  </span>
                )}
              </div>

              <div
                className={`max-w-[88%] rounded-2xl p-3.5 text-xs leading-relaxed ${
                  m.role === 'user'
                    ? 'bg-gradient-to-r from-cyan-600 to-indigo-600 text-white rounded-tr-none shadow-md'
                    : 'bg-slate-800/90 text-slate-100 border border-slate-700/80 rounded-tl-none shadow-sm'
                }`}
              >
                <div className="whitespace-pre-wrap font-sans selection:bg-cyan-500/30">
                  {m.text}
                </div>

                {detectedJson && (
                  <div className="mt-3 pt-3 border-t border-slate-700/60 flex items-center justify-between gap-2">
                    <span className="text-[10px] text-cyan-400 font-mono font-semibold flex items-center gap-1">
                      <Check className="w-3 h-3" />
                      JSON gotowy do pobrania
                    </span>
                    <div className="flex items-center gap-1.5">
                      <button
                        onClick={() => handleCopy(detectedJson, m.id)}
                        className="px-2 py-1 bg-slate-700 hover:bg-slate-600 text-slate-200 rounded text-[10px] font-mono flex items-center gap-1 transition cursor-pointer"
                      >
                        {copiedId === m.id ? <Check className="w-3 h-3 text-emerald-400" /> : <Copy className="w-3 h-3" />}
                        <span>{copiedId === m.id ? 'Skopiowano' : 'Kopiuj'}</span>
                      </button>
                      <button
                        onClick={() => handleDownloadJson(detectedJson)}
                        className="px-2 py-1 bg-cyan-700 hover:bg-cyan-600 text-white rounded text-[10px] font-mono flex items-center gap-1 transition cursor-pointer"
                      >
                        <Download className="w-3 h-3" />
                        <span>Pobierz .json</span>
                      </button>
                    </div>
                  </div>
                )}
              </div>
            </div>
          );
        })}

        {isLoading && (
          <div className="flex items-center gap-2 text-cyan-400 font-mono text-xs p-2 bg-slate-950/40 rounded-xl border border-slate-800/60 w-fit">
            <Sparkles className="w-4 h-4 animate-spin text-cyan-300" />
            <span>Obliczanie trajektorii i kwantyzacja E8...</span>
          </div>
        )}
        <div ref={messagesEndRef} />
      </div>

      {/* Quick Prompts Bar */}
      <div className="px-3 py-2 bg-slate-950/80 border-t border-slate-800/70 overflow-x-auto shrink-0">
        <div className="flex items-center gap-1.5 whitespace-nowrap">
          <span className="text-[10px] font-mono text-slate-400 shrink-0 font-bold uppercase">Szybkie pytania:</span>
          {quickPrompts.map((q, idx) => (
            <button
              key={idx}
              onClick={() => handleSend(q)}
              disabled={isLoading}
              className="text-[10px] font-mono px-2 py-1 bg-slate-800/90 hover:bg-slate-750 border border-slate-700 text-slate-300 hover:text-white rounded-lg transition-colors cursor-pointer shrink-0 disabled:opacity-40"
            >
              {idx === 0 ? "★ amino_acid_warehouse.json (GLY, ALA, VAL)" : q.slice(0, 36) + "..."}
            </button>
          ))}
        </div>
      </div>

      {/* Input field */}
      <div className="p-3 bg-slate-950 border-t border-slate-800 shrink-0">
        <form
          onSubmit={(e) => {
            e.preventDefault();
            handleSend();
          }}
          className="flex items-center gap-2"
        >
          <input
            type="text"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            placeholder="Zapytaj architekta SNAPLOCK (np. 'wygeneruj warehouse dla GLY, ALA, VAL')..."
            className="flex-1 bg-slate-900 border border-slate-700/80 focus:border-cyan-500 rounded-xl px-3.5 py-2 text-xs text-slate-100 placeholder-slate-500 focus:outline-none transition-colors"
            disabled={isLoading}
          />
          <button
            type="submit"
            disabled={!input.trim() || isLoading}
            className="px-4 py-2 bg-gradient-to-r from-cyan-600 to-indigo-600 hover:from-cyan-500 hover:to-indigo-500 disabled:opacity-40 text-white rounded-xl text-xs font-semibold flex items-center gap-1.5 shadow-md transition-all cursor-pointer shrink-0"
          >
            <Send className="w-3.5 h-3.5" />
            <span>Wyślij</span>
          </button>
          <button
            type="button"
            onClick={() => setMessages([messages[0]])}
            className="p-2 text-slate-400 hover:text-rose-400 hover:bg-slate-900 rounded-xl transition-colors cursor-pointer shrink-0"
            title="Wyczyść historię"
          >
            <Trash2 className="w-4 h-4" />
          </button>
        </form>
      </div>
    </div>
  );
};
