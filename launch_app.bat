@echo off
setlocal
set "ROOT=%~dp0"
cd /d "%ROOT%"

set "VENV_PY=%ROOT%.venv\Scripts\python.exe"
if exist "%VENV_PY%" (
	"%VENV_PY%" -m streamlit run examples\gui_app.py
) else (
	python -m streamlit run examples\gui_app.py
)
