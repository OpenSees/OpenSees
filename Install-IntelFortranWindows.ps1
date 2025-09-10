$Url = "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/2bff9d18-ac2b-40b1-8167-00156f466b0e/intel-fortran-essentials-2025.0.0.301_offline.exe"
$Hash = "E1DF578C069E567D14E0EA124DA3AC921877CA19B46D1251F6D18A3D25F0CE615F496FAE75095DEC83251A9D2E0459E1"
$Components = "intel.oneapi.win.ifort-compiler:intel.oneapi.win.mkl.devel:intel.oneapi.win.mpi.devel"

$FullVersion = $Url.Split("-")[-1].Split("_")[0]
$Version = $FullVersion.Split(".")[0..1] -join "."

function Test-IntelInstallation {
    param (
        [string]$Version
    )
    $PathsToCheck = @(
        "C:\Program Files (x86)\Intel\oneAPI\$Version\bin\ifx.exe",
        "C:\Program Files (x86)\Intel\oneAPI\$Version\bin\mpiexec.exe",
        "C:\Program Files (x86)\Intel\oneAPI\$Version\bin\mkl_intel_thread.2.dll"
    )

    $Missing = @()
    foreach ($Path in $PathsToCheck) {
        if (-not (Test-Path -Path $Path)) {
            $Missing += $Path
        }
    }

    if ($Missing.Count -gt 0) {
        Write-Host "The following files are missing:" -ForegroundColor Yellow
        foreach ($File in $Missing) {
            Write-Host "  $File" -ForegroundColor Yellow
        }
        return $false
    } else {
        Write-Host "All checked components are installed." -ForegroundColor Green
        return $true
    }
}

if (Test-IntelInstallation -Version $Version) {
    Write-Host "Intel Fortran Essentials $Version is already installed." -ForegroundColor Green
    exit 0
}

$DownloadDir = $PSScriptRoot + "/downloads"
$Filename = $Url.Split("/")[-1]
if (!(Test-Path -Path $DownloadDir)) {
    New-Item -ItemType Directory -Path $DownloadDir | Out-Null
}

$OutputPath = Join-Path -Path $DownloadDir -ChildPath $Filename

if ((Test-Path -Path $OutputPath) -and (Get-FileHash -Path $OutputPath -Algorithm SHA384).Hash -eq $Hash) {
    Write-Host "File already downloaded and verified: $OutputPath"
} else {
    Invoke-WebRequest -Uri $Url -OutFile $OutputPath
    if ((Get-FileHash -Path $OutputPath -Algorithm SHA384).Hash -ne $Hash) {
        throw "Downloaded file hash does not match expected hash."
    }
}

$ExtractedPath = $OutputPath + "_extracted"
if (Test-Path -Path $ExtractedPath) {
    Remove-Item -Recurse -Force -Path $ExtractedPath
}
$LogPath = $OutputPath + "_extract.log"

$ExtractArgs = @("-s", "-x", "-f", $ExtractedPath, "--log", $LogPath)
$Process = Start-Process -FilePath $OutputPath -ArgumentList $ExtractArgs -Wait -PassThru
if ($Process.ExitCode -ne 0) {
    throw "Extraction failed with exit code $($Process.ExitCode)"
}

$BootstrapperPath = Join-Path -Path $ExtractedPath -ChildPath "bootstrapper.exe"

$LogDir = $OutputPath + "_install_logs"
if (!(Test-Path -Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$Process = Start-Process -FilePath $BootstrapperPath -ArgumentList "-s", "--action", "install", "--eula=accept", "--log-dir=$LogDir", "--components=$Components" -Wait -PassThru
if ($Process.ExitCode -ne 0) {
    throw "Installation failed with exit code $($Process.ExitCode)"
}

if (-not (Test-IntelInstallation -Version $Version)) {
    throw "Installation verification failed."
}
