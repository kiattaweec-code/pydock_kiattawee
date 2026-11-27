; Inno Setup script for pydock.
;
;   iscc /DMyAppVersion=0.2.1 packaging\windows_installer.iss
;
; Installs per-user by default, so no administrator prompt is needed. The app
; writes its jobs and kernel cache to %LOCALAPPDATA%\pydock, which this
; removes on uninstall only if the user asks.

#ifndef MyAppVersion
  #define MyAppVersion "0.2.1"
#endif

#define MyAppName "pydock"
#define MyAppExeName "pydock.exe"

[Setup]
AppId={{8E4C2A91-6D3F-4B57-9E1A-7C2D5F8B4A33}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher=pydock
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
OutputDir=..\dist
OutputBaseFilename=pydock-{#MyAppVersion}-windows-x86_64-setup
SetupIconFile=icon.ico
UninstallDisplayIcon={app}\{#MyAppExeName}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
; Per-user by default: no UAC prompt, and the app needs nothing privileged.
PrivilegesRequiredOverridesAllowed=dialog
PrivilegesRequired=lowest

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; \
  GroupDescription: "Additional shortcuts:"

[Files]
Source: "..\dist\pydock\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; \
  Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; \
  Description: "Start {#MyAppName}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Leave the user's docking results alone unless they say otherwise; this only
; removes the compiled-kernel cache, which is regenerated on next launch.
Type: filesandordirs; Name: "{localappdata}\pydock\cache"

[Code]
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  DataDir: String;
begin
  if CurUninstallStep = usPostUninstall then
  begin
    DataDir := ExpandConstant('{localappdata}\pydock');
    if DirExists(DataDir) then
      if MsgBox('Also delete your saved docking jobs in' + #13#10 + DataDir + '?',
                mbConfirmation, MB_YESNO) = IDYES then
        DelTree(DataDir, True, True, True);
  end;
end;
