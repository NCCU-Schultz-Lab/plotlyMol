Set shell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
root = fso.GetParentFolderName(WScript.ScriptFullName)
batPath = Chr(34) & root & "\stop_app.bat" & Chr(34)
shell.Run batPath, 0
