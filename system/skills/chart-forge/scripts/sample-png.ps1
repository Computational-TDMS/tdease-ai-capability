param([string]$Png, [double[]]$Stops = @(0.0, 0.12, 0.3, 0.42, 0.55, 0.75, 1.0))
Add-Type -AssemblyName System.Drawing
$img = [System.Drawing.Bitmap]::FromFile($Png)
Write-Host "size: $($img.Width)x$($img.Height)"
for ($b = 0; $b -lt $Stops.Count - 1; $b++) {
  $y0 = [int]($img.Height * $Stops[$b]); $y1 = [int]($img.Height * $Stops[$b + 1])
  if ($y1 -le $y0) { continue }
  $colored = 0; $dark = 0; $mid = 0
  for ($y = $y0; $y -lt $y1; $y += 3) {
    for ($x = 8; $x -lt $img.Width - 8; $x += 3) {
      $p = $img.GetPixel($x, $y)
      $mx = [Math]::Max($p.R, [Math]::Max($p.G, $p.B))
      $mn = [Math]::Min($p.R, [Math]::Min($p.G, $p.B))
      if ($mx - $mn -gt 40) { $colored++ }
      if ($mx -lt 120) { $dark++ }
      if ($mn -ge 110 -and $mx -le 210) { $mid++ }
    }
  }
  Write-Host ("band {0,3}-{1,3}%: colored {2,5}  dark {3,5}  mid {4,5}" -f [int]($Stops[$b]*100), [int]($Stops[$b+1]*100), $colored, $dark, $mid)
}
$img.Dispose()
