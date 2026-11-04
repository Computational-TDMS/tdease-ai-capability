param([string]$Png, [int]$Cols = 110, [int]$Rows = 56)
Add-Type -AssemblyName System.Drawing
$img = [System.Drawing.Bitmap]::FromFile($Png)
$chars = " .:-=+*#%@"
for ($r = 0; $r -lt $Rows; $r++) {
  $line = ""
  for ($c = 0; $c -lt $Cols; $c++) {
    $x0 = [int]($img.Width * $c / $Cols); $x1 = [int]($img.Width * ($c + 1) / $Cols)
    $y0 = [int]($img.Height * $r / $Rows); $y1 = [int]($img.Height * ($r + 1) / $Rows)
    $dark = 0; $total = 0
    for ($y = $y0; $y -lt $y1; $y += 2) {
      for ($x = $x0; $x -lt $x1; $x += 2) {
        $p = $img.GetPixel($x, $y)
        $mx = [Math]::Max($p.R, [Math]::Max($p.G, $p.B))
        if ($mx -lt 200) { $dark++ }
        $total++
      }
    }
    $d = if ($total -gt 0) { $dark / $total } else { 0 }
    $idx = [int]([Math]::Min(9, $d * 9 * 2.2))
    $line += $chars[$idx]
  }
  Write-Host $line
}
$img.Dispose()
