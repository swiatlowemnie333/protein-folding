import express from 'express';
import fs from 'fs';
import { createServer as createViteServer } from 'vite';
import path from 'path';
import { fileURLToPath } from 'url';
import dotenv from 'dotenv';
import { GoogleGenAI } from '@google/genai';
import { 
  SNAPLOCK_22_DATABASE, 
  generateAminoAcidWarehouseJson, 
  evaluateSnaplockJunction,
  runSnaplockCalibrationFromPdb
} from './src/snaplockEngine.ts';

dotenv.config();

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const SNAPLOCK_SYSTEM_INSTRUCTION = `Jesteś głównym inżynierem biofizyki obliczeniowej i architektem deterministycznego silnika zwijania białek o nazwie SNAPLOCK. Twoim celem jest rezygnacja z losowego przeszukiwania ciągłych kątów torsyjnych na rzecz dyskretnego mapowania geometrycznego opartego na sferycznej siatce 248 wektorów (algebrze E8).

OFICJALNY MODUŁ WALIDACYJNY I BENCHMARK ALPHAFOLD (snaplock_alphafold_benchmark.py):
Każda proponowana struktura przechodzi kwantyzację przez siatkę 248 wektorów:
1. Generowanie siatki 248 wektorów E8 (kardynalne: #1, #248, #100, #200, #50, #150 + spirala Fibonacciego).
2. Wektory wiązań C_alfa_i -> C_alfa_{i+1} rzutowane są na najbliższą łatkę: Patch_best = argmax (v_ref . e_i).
3. Obliczanie średniego błędu kątowego siatki (mean_angular_patch_loss_deg, typowo ok. 5.0°).
4. Wykrywanie klamer kowalencyjnych S-S Cys-Cys (SG-SG <= 2.5 Å) z najwyższym priorytetem.
5. Zgodność z 5 Ścieżkami SNAPLOCK (szkielet #248->#1, wykluczenie kolizji VdW, elektrostatyka pH 7.4, pakowanie rdzenia, mostki kowalencyjne).

Mechanika SNAPLOCK (Żelazne Reguły):
1. Każdy z 21 aminokwasów jest wpisany w sztywną sferę, na której powierzchni zdefiniowano 248 dyskretnych "łatek" kontaktowych.
2. Układ lokalny każdego aminokwasu (środek w punkcie geometrycznym szkieletu (0,0,0), oparty na prawoskrętnym układzie kartezjańskim 3D):
   - Góra (+Y): Łatka #1 (Wejście: atom N, port wiązania peptydowego) -> wektor [0, 1, 0].
   - Dół (-Y): Łatka #248 (Wyjście: atom C, tlen karbonylowy) -> wektor [0, -1, 0].
   - Przód (+Z): Łatka #50 -> wektor [0, 0, 1].
   - Prawo (+X): Łatka #100 -> wektor [1, 0, 0].
   - Tył (-Z): Łatka #150 -> wektor [0, 0, -1].
   - Lewo (-X): Łatka #200 -> wektor [-1, 0, 0].
3. Dystans między środkami układów (0,0,0) sąsiadujących aminokwasów służy do natychmiastowej detekcji kolizji i "zgniecenia" na podstawie ich fizycznej objętości Van der Waalsa.
4. Aminokwasy łączą się ze sobą wyłącznie łatkami #248 (poprzednik) i #1 (następnik). Jedynym stopniem swobody jest dyskretna rotacja wokół osi styku w zakresie od -180° do +180°, skacząca po sąsiadujących łatkach.
5. Ewaluacja złączy odbywa się na podstawie 5 rygorystycznych ścieżek (Tracks):
   - Track 1: weryfikacja szkieletu (backbone porty #248 -> #1),
   - Track 2: wykluczenie kolizji Van der Waalsa (wykluczenie zgniecenia sterycznego),
   - Track 3: optymalizacja ładunków w pH 7.4 (mostki solne i odpychanie jednoimienne),
   - Track 4: zamknięcie rdzenia hydrofobowego (hydropatia Kyte-Doolittle i desolwatacja),
   - Track 5: klamry kowalencyjne (mostki dwusiarczkowe S-S CYS-CYS / seleniaczkowe Se-Se/Se-S SEC-SEC/SEC-CYS nadpisujące słabsze oddziaływania z najwyższym priorytetem).

Baza Danych SNAPLOCK (22 Parametry Fizykochemiczne):
Dysponujesz dokładną matrycą 22 parametrów dla 21 aminokwasów (GLY, ALA, VAL, LEU, ILE, PRO, PHE, TYR, TRP, SER, THR, CYS, MET, ASN, GLN, ASP, GLU, LYS, ARG, HIS, SEC).
Gdy użytkownik pyta o wygenerowanie pliku amino_acid_warehouse.json lub walidację AlphaFold, dostarcz profesjonalną, sformatowaną strukturę JSON oraz dokładne wyliczenia zgodne z oficjalnym protokołem.`;

async function startServer() {
  const app = express();
  const PORT = Number(process.env.PORT) || 3000;

  app.use(express.json({ limit: '10mb' }));

  // Initialize server-side Gemini client
  const ai = new GoogleGenAI({
    apiKey: process.env.GEMINI_API_KEY || '',
    httpOptions: {
      headers: {
        'User-Agent': 'aistudio-build',
      }
    }
  });

  // API Route: Multi-turn Gemini chat for SNAPLOCK
  app.post('/api/gemini/chat', async (req, res) => {
    try {
      const { messages, modelName = 'gemini-3.8-flash' } = req.body;
      if (!messages || !Array.isArray(messages) || messages.length === 0) {
        return res.status(400).json({ error: 'Brak wiadomości w historii czatu.' });
      }

      // Check if API key is present
      if (!process.env.GEMINI_API_KEY) {
        // Fallback deterministic engineer response if key not configured
        const lastMsg = messages[messages.length - 1].text || '';
        let fallbackResponse = '';
        if (lastMsg.toLowerCase().includes('amino_acid_warehouse.json') || lastMsg.toLowerCase().includes('warehouse')) {
          const warehouse = generateAminoAcidWarehouseJson(['GLY', 'ALA', 'VAL']);
          fallbackResponse = `Raport Głównego Inżyniera SNAPLOCK:
Wygenerowano strukturę magazynu \`amino_acid_warehouse.json\` dla pierwszych trzech aminokwasów (GLY, ALA, VAL) z pełnym odwzorowaniem 22 parametrów fizykochemicznych na dyskretną siatkę 248 wektorów algebry E8.

\`\`\`json
${JSON.stringify(warehouse, null, 2)}
\`\`\`

**Kluczowe właściwości geometryczne E8:**
- **GLY (#1)**: Objętość VdW 48 Å³, brak rotamerów, symetryczne wejście/wyjście.
- **ALA (#2)**: Objętość VdW 67 Å³, preferencja helisy alfa, promień VdW = 2.52 Å.
- **VAL (#3)**: Objętość VdW 105 Å³, silna hydropatia Kyte +4.2, domknięcie rdzenia hydrofobowego (Track 4).

Wszystkie połączenia zachowują sztywne złącze #248 (poprzednik) -> #1 (następnik) ze skokową rotacją wokół osi wiązania peptydowego.`;
        } else {
          fallbackResponse = `Jako Główny Inżynier Biofizyki Obliczeniowej SNAPLOCK melduję:
Silnik operuje deterministycznie na dyskretnej siatce 248 wektorów algebry E8. Złącza aminokwasowe są weryfikowane według 5 ścieżek ewaluacyjnych (Track 1-5). Wprowadź sekwencję lub zapytaj o parametry magazynu 21 aminokwasów.`;
        }
        return res.json({ text: fallbackResponse });
      }

      // Map format to Gemini API contents
      // System instruction passed via config.systemInstruction
      const contents = messages.map((m: { role: string; text: string }) => ({
        role: m.role === 'model' || m.role === 'assistant' ? 'model' : 'user',
        parts: [{ text: m.text }]
      }));

      // Selected valid model from supported list
      const allowedModels = ['gemini-4-argon', 'gemini-3.8-flash', 'gemini-3.1-flash-lite'];
      const targetModel = allowedModels.includes(modelName) ? modelName : 'gemini-3.8-flash';

      const response = await ai.models.generateContent({
        model: targetModel,
        contents,
        config: {
          systemInstruction: SNAPLOCK_SYSTEM_INSTRUCTION,
          temperature: 0.2, // deterministic, precise scientific reasoning
        }
      });

      const responseText = response.text || 'Brak odpowiedzi od silnika Gemini.';
      return res.json({ text: responseText, model: targetModel });
    } catch (err: any) {
      console.error('Error in /api/gemini/chat:', err);
      // Deterministic fallback if API error occurs
      const warehouse = generateAminoAcidWarehouseJson(['GLY', 'ALA', 'VAL']);
      return res.json({
        text: `Wygenerowano dane w trybie bezpośrednim silnika SNAPLOCK (E8):\n\n\`\`\`json\n${JSON.stringify(warehouse, null, 2)}\n\`\`\``,
        warning: err.message || 'Wystąpił błąd API Gemini, użyto silnika lokalnego.'
      });
    }
  });

  // API Route: Direct generation of amino_acid_warehouse.json
  app.post('/api/gemini/warehouse', (req, res) => {
    try {
      const { codes = ['GLY', 'ALA', 'VAL'] } = req.body;
      const data = generateAminoAcidWarehouseJson(codes);
      res.json(data);
    } catch (e: any) {
      res.status(500).json({ error: e.message });
    }
  });

  // API Route: Evaluate 5 Tracks for arbitrary pair
  app.post('/api/snaplock/evaluate-junction', (req, res) => {
    try {
      const { predCode, succCode, torsion = 0 } = req.body;
      const pred = SNAPLOCK_22_DATABASE.find(a => a.code3 === predCode) || SNAPLOCK_22_DATABASE[0];
      const succ = SNAPLOCK_22_DATABASE.find(a => a.code3 === succCode) || SNAPLOCK_22_DATABASE[1];
      const result = evaluateSnaplockJunction(pred, succ, 248, 1, torsion);
      res.json(result);
    } catch (e: any) {
      res.status(500).json({ error: e.message });
    }
  });

  // API Route: Proxy AlphaFold DB prediction metadata (Science Skill requirement)
  app.get('/api/alphafold/prediction/:uniprotId', async (req, res) => {
    try {
      const { uniprotId } = req.params;
      const response = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);
      if (!response.ok) {
        return res.status(response.status).json({ error: `AlphaFold API error: ${response.statusText}` });
      }
      const data = await response.json();
      res.json(data);
    } catch (err: any) {
      res.status(500).json({ error: err.message });
    }
  });

  // API Route: SNAPLOCK vs AlphaFold Benchmark & Calibration Engine
  app.post('/api/snaplock/benchmark', async (req, res) => {
    try {
      const { uniprotId = 'P83347', pdbContent } = req.body;
      let targetPdb = pdbContent;

      if (!targetPdb && uniprotId) {
        const localPath = path.resolve(__dirname, `${uniprotId}_alphafold.pdb`);
        if (fs.existsSync(localPath)) {
          targetPdb = fs.readFileSync(localPath, 'utf-8');
        } else {
          // Fetch from AlphaFold DB EMBL-EBI
          try {
            const metaResp = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);
            if (metaResp.ok) {
              const metaData = await metaResp.json();
              if (metaData && metaData[0] && metaData[0].pdbUrl) {
                const pdbResp = await fetch(metaData[0].pdbUrl);
                if (pdbResp.ok) {
                  targetPdb = await pdbResp.text();
                }
              }
            }
          } catch (fetchErr) {
            console.warn(`Could not fetch online PDB for ${uniprotId}, checking local files...`);
          }
        }
      }

      if (!targetPdb) {
        return res.status(404).json({ error: `Nie odnaleziono modelu AlphaFold dla ID: ${uniprotId}` });
      }

      const report = runSnaplockCalibrationFromPdb(targetPdb, uniprotId);
      res.json(report);
    } catch (err: any) {
      console.error('Error in /api/snaplock/benchmark:', err);
      res.status(500).json({ error: err.message });
    }
  });

  // ESMFold (Meta) - prawdziwa predykcja struktury z sekwencji
  app.post('/api/esmfold', async (req, res) => {
    const seq = String(req.body?.sequence || '').toUpperCase().replace(/[^A-Z]/g, '');
    if (seq.length < 10 || seq.length > 400 || /[^ACDEFGHIKLMNPQRSTVWY]/.test(seq)) {
      return res.status(400).json({ error: 'Sekwencja: 10-400 standardowych aminokwasow' });
    }
    try {
      const r = await fetch('https://api.esmatlas.com/foldSequence/v1/pdb/', {
        method: 'POST',
        headers: { 'Content-Type': 'text/plain' },
        body: seq,
        signal: AbortSignal.timeout(120000)
      });
      const pdb = await r.text();
      if (!r.ok || !pdb.includes('ATOM')) {
        return res.status(502).json({ error: `ESMFold HTTP ${r.status}` });
      }
      res.json({ pdb });
    } catch (err: any) {
      res.status(502).json({ error: `ESMFold niedostepny: ${err.message}` });
    }
  });

  // Health check endpoints for Cloud Run
  app.get(['/healthz', '/api/health'], (_req, res) => {
    res.status(200).json({ status: 'ok', engine: 'SNAPLOCK E8' });
  });

  // Vite middleware in dev or static files in production
  const distHtmlPath = path.resolve(__dirname, 'dist', 'index.html');
  const isDev = process.env.npm_lifecycle_event === 'dev' || (!fs.existsSync(distHtmlPath) && process.env.NODE_ENV !== 'production');
  const isProduction = !isDev;

  if (isProduction) {
    console.log('[SNAPLOCK Engine] Serwowanie plików produkcyjnych z dist/');
    app.use(express.static(path.resolve(__dirname, 'dist')));
    app.get('*', (_req, res) => {
      res.sendFile(distHtmlPath);
    });
  } else {
    console.log('[SNAPLOCK Engine] Uruchamianie Vite middleware dla trybu deweloperskiego');
    const vite = await createViteServer({
      server: {
        middlewareMode: true,
        host: '0.0.0.0',
        hmr: false
      },
      appType: 'spa',
    });
    app.use(vite.middlewares);
  }

  app.listen(PORT, '0.0.0.0', () => {
    console.log(`[SNAPLOCK Engine] Serwer uruchomiony na http://0.0.0.0:${PORT}`);
  });
}

startServer();
