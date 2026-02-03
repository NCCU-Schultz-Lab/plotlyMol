@echo off
setlocal
set "ROOT=%~dp0"
cd /d "%ROOT%"

powershell -NoProfile -Command "Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -like '*streamlit*run*examples\\gui_app.py*' } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force }"
