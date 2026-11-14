# download_pdb.ps1 - pobiera niepowtarzalne bialka z RCSB (X-ray <=2.0 A, 50-300 aa, klastry 30% identycznosci)
# i zapisuje TYLKO atomy N,CA,C pierwszego lancucha (male pliki). Na koncu pakuje do bialka_backbone.zip.
# Uruchom: PowerShell -ExecutionPolicy Bypass -File download_pdb.ps1 [-Count 1500] [-OutDir C:\bialka]
# Nic nie instaluje, nie wymaga hasel. Mozna przerwac (Ctrl+C) i uruchomic ponownie - pomija pobrane.
param([int]$Count = 1500, [string]$OutDir = "$PSScriptRoot\bialka_pdb")
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

Write-Host "1/3 Pobieram liste bialek z RCSB (Count=$Count) ..."
$q = @{
  query = @{ type='group'; logical_operator='and'; nodes=@(
    @{ type='terminal'; service='text'; parameters=@{ attribute='exptl.method'; operator='exact_match'; value='X-RAY DIFFRACTION' } },
    @{ type='terminal'; service='text'; parameters=@{ attribute='rcsb_entry_info.resolution_combined'; operator='less_or_equal'; value=2.0 } },
    @{ type='terminal'; service='text'; parameters=@{ attribute='rcsb_entry_info.polymer_entity_count_protein'; operator='equals'; value=1 } },
    @{ type='terminal'; service='text'; parameters=@{ attribute='rcsb_entry_info.deposited_polymer_monomer_count'; operator='range'; value=@{ from=50; to=300; include_lower=$true; include_upper=$true } } }
  ) }
  return_type = 'polymer_entity'
  request_options = @{
    group_by = @{ aggregation_method='sequence_identity'; similarity_cutoff=30; ranking_criteria_type=@{ sort_by='rcsb_entry_info.resolution_combined'; direction='asc' } }
    group_by_return_type = 'representatives'
    paginate = @{ start=0; rows=$Count }
    results_verbosity = 'minimal'
  }
}
$json = $q | ConvertTo-Json -Depth 12
try {
  $r = Invoke-RestMethod -Method Post -Uri 'https://search.rcsb.org/rcsbsearch/v2/query' -ContentType 'application/json' -Body $json
} catch {
  Write-Host "Blad zapytania do RCSB: $($_.Exception.Message)"
  Write-Host "Awaryjnie: wpisz ID (np. 1UBQ) po jednym w linii do pliku ids.txt obok skryptu i uruchom ponownie."
  if (-not (Test-Path "$PSScriptRoot\ids.txt")) { exit 1 }
  $r = $null
}
if ($r) { $ids = $r.result_set | ForEach-Object { ($_.identifier -split '_')[0] } | Select-Object -Unique }
else    { $ids = Get-Content "$PSScriptRoot\ids.txt" | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim().ToUpper() } }
Write-Host "Znaleziono $($ids.Count) struktur."

$std = 'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
$wc = New-Object System.Net.WebClient
$ok = 0; $skip = 0; $fail = 0; $i = 0
Write-Host "2/3 Pobieram struktury ..."
foreach ($id in $ids) {
  $i++
  $out = Join-Path $OutDir "$id.pdb"
  if (Test-Path $out) { $ok++; continue }
  try {
    $txt = $wc.DownloadString("https://files.rcsb.org/download/$id.pdb")
    $chain = $null; $keep = New-Object System.Collections.Generic.List[string]; $res = @{}
    foreach ($line in ($txt -split "`n")) {
      if ($line.StartsWith('ENDMDL')) { break }
      if (-not $line.StartsWith('ATOM')) { continue }
      if ($line.Length -lt 54) { continue }
      if ($line[16] -ne ' ' -and $line[16] -ne 'A') { continue }
      $ch = $line[21]; if ($null -eq $chain) { $chain = $ch }; if ($ch -ne $chain) { continue }
      $rn = $line.Substring(17,3); if ($std -notcontains $rn) { continue }
      $an = $line.Substring(12,4).Trim(); if ($an -ne 'N' -and $an -ne 'CA' -and $an -ne 'C') { continue }
      $keep.Add($line.Substring(0,54)); $res[$line.Substring(22,5)] = 1
    }
    if ($res.Count -ge 40) { [IO.File]::WriteAllLines($out, $keep); $ok++ } else { $skip++ }
  } catch { $fail++ }
  if ($i % 50 -eq 0) { Write-Host "  $i/$($ids.Count)  ok=$ok pominiete=$skip bledy=$fail" }
}
Write-Host "3/3 Pakuje ..."
$zip = Join-Path $PSScriptRoot 'bialka_backbone.zip'
if (Test-Path $zip) { Remove-Item $zip }
Compress-Archive -Path "$OutDir\*.pdb" -DestinationPath $zip
$mb = [math]::Round((Get-Item $zip).Length / 1MB, 1)
Write-Host "GOTOWE: $ok struktur (pominiete=$skip, bledy=$fail). Zip: $zip ($mb MB)"
Write-Host "Wrzuc ten zip do rozmowy z Claude."
