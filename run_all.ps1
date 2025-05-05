# Step 1: Clean
Write-Host "Cleaning"
make clean

# Step 2: Compile
Write-Host "`nCompiling"
make

# Step 3: Delete old output
if (Test-Path "tst.out") {
    Remove-Item "tst.out"
    Write-Host "`nDeleted old tst.out"
}

# Step 4: Launch Fortran (teach.exe) in background and capture its process
Write-Host "`nRunning teach.exe < Input_fine in background..."
$fortranProcess = Start-Process powershell -NoNewWindow -PassThru -ArgumentList "-Command", "Get-Content Input_fine | .\teach.exe"

# Step 5: Wait briefly for file to start writing
Start-Sleep -Seconds 1

# Step 6: Run Python plot in foreground (blocks until you close it)
Write-Host "`nLaunching Live_Residuals.py..."
python Live_Residuals.py

# Step 7: After plot closes, check Fortran status
if ($fortranProcess.HasExited) {
    Write-Host "`nteach.exe already finished."
} else {
    Write-Host "teach.exe still running... waiting for it to finish."
    Wait-Process -Id $fortranProcess.Id
    Write-Host "teach.exe has now finished."
}
