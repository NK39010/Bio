param(
    [string]$ImageName = "orivfinder-ready",
    [string]$ImageTag = "latest",
    [string]$InputDir = "data/input",
    [string]$OutputDir = "data/output",
    [string[]]$Patterns = @("*.fasta", "*.fa", "*.fna"),
    [switch]$LoadImage,
    [string]$ImageArchive = "orivfinder-ready.tar.gz",
    [switch]$Overwrite,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

function Resolve-ProjectPath {
    param([string]$PathText)
    if ([System.IO.Path]::IsPathRooted($PathText)) {
        return (Resolve-Path -LiteralPath $PathText).Path
    }
    return (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot $PathText)).Path
}

function Ensure-Directory {
    param([string]$PathText)
    if (-not (Test-Path -LiteralPath $PathText)) {
        New-Item -ItemType Directory -Path $PathText -Force | Out-Null
    }
}

function To-DockerPath {
    param(
        [string]$FullPath,
        [string]$MountedRoot,
        [string]$DockerRoot = "/app/data"
    )

    $relative = [System.IO.Path]::GetRelativePath($MountedRoot, $FullPath)
    $relative = $relative -replace "\\", "/"
    return "$DockerRoot/$relative"
}

$docker = Get-Command docker -ErrorAction SilentlyContinue
if (-not $docker) {
    throw "Docker was not found on PATH. Please install/start Docker Desktop first."
}

$projectRoot = $PSScriptRoot
$dataRoot = Join-Path $projectRoot "data"
if ([System.IO.Path]::IsPathRooted($InputDir)) {
    $inputRoot = $InputDir
}
else {
    $inputRoot = Join-Path $projectRoot $InputDir
}
if ([System.IO.Path]::IsPathRooted($OutputDir)) {
    $outputRoot = $OutputDir
}
else {
    $outputRoot = Join-Path $projectRoot $OutputDir
}

Ensure-Directory $dataRoot
Ensure-Directory $inputRoot
Ensure-Directory $outputRoot

if ($LoadImage) {
    $archivePath = Join-Path $projectRoot $ImageArchive
    if (-not (Test-Path -LiteralPath $archivePath)) {
        throw "Image archive not found: $archivePath"
    }

    Write-Host "Loading Docker image from $archivePath"
    if ($archivePath.EndsWith(".gz")) {
        if (Get-Command gzip -ErrorAction SilentlyContinue) {
            & gzip -dc $archivePath | docker load
        }
        else {
            throw "gzip was not found. Extract $ImageArchive to .tar first, then run: docker load -i <file.tar>"
        }
    }
    else {
        docker load -i $archivePath
    }
}

$imageRef = "${ImageName}:${ImageTag}"
Write-Host "Checking Docker image: $imageRef"
$imageFound = docker images --format "{{.Repository}}:{{.Tag}}" | Where-Object { $_ -eq $imageRef }
if (-not $imageFound) {
    throw "Docker image '$imageRef' was not found. Run with -LoadImage first, or load it manually."
}

$inputFiles = New-Object System.Collections.Generic.List[System.IO.FileInfo]
foreach ($pattern in $Patterns) {
    Get-ChildItem -LiteralPath $inputRoot -File -Filter $pattern | ForEach-Object {
        $inputFiles.Add($_)
    }
}

$inputFiles = $inputFiles | Sort-Object FullName -Unique
if (-not $inputFiles -or $inputFiles.Count -eq 0) {
    Write-Host "No input FASTA files found in: $inputRoot"
    Write-Host "Put .fasta, .fa, or .fna files under data/input, then run this script again."
    exit 0
}

$summaryPath = Join-Path $outputRoot "batch_summary.tsv"
"input_file`toutput_dir`tstatus`tstarted_at`tfinished_at" | Set-Content -LiteralPath $summaryPath -Encoding UTF8

foreach ($file in $inputFiles) {
    $name = [System.IO.Path]::GetFileNameWithoutExtension($file.Name)
    $sampleOutput = Join-Path $outputRoot $name
    $doneMarker = Join-Path $sampleOutput ".done"

    if ((Test-Path -LiteralPath $doneMarker) -and -not $Overwrite) {
        Write-Host "Skipping completed input: $($file.Name)"
        "$($file.FullName)`t$sampleOutput`tskipped_existing`t`t" | Add-Content -LiteralPath $summaryPath -Encoding UTF8
        continue
    }

    if ((Test-Path -LiteralPath $sampleOutput) -and $Overwrite) {
        Remove-Item -LiteralPath $sampleOutput -Recurse -Force
    }
    Ensure-Directory $sampleOutput

    $dockerInput = "/app/input/$($file.Name)"
    $dockerOutput = "/app/output/$name"
    $started = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

    $dockerArgs = @(
        "run",
        "--rm",
        "-v",
        "${inputRoot}:/app/input:ro",
        "-v",
        "${outputRoot}:/app/output",
        $imageRef,
        "python",
        "oriVfinder.py",
        "--fasta",
        $dockerInput,
        "--output_dir",
        $dockerOutput
    )

    Write-Host ""
    Write-Host "Input:  $($file.FullName)"
    Write-Host "Output: $sampleOutput"
    Write-Host "Docker: docker $($dockerArgs -join ' ')"

    if ($DryRun) {
        "$($file.FullName)`t$sampleOutput`tdry_run`t$started`t" | Add-Content -LiteralPath $summaryPath -Encoding UTF8
        continue
    }

    try {
        & docker @dockerArgs
        if ($LASTEXITCODE -ne 0) {
            throw "docker exited with code $LASTEXITCODE"
        }
        New-Item -ItemType File -Path $doneMarker -Force | Out-Null
        $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        "$($file.FullName)`t$sampleOutput`tsuccess`t$started`t$finished" | Add-Content -LiteralPath $summaryPath -Encoding UTF8
    }
    catch {
        $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        "$($file.FullName)`t$sampleOutput`tfailed: $($_.Exception.Message)`t$started`t$finished" | Add-Content -LiteralPath $summaryPath -Encoding UTF8
        Write-Error "Failed on $($file.Name): $($_.Exception.Message)"
    }
}

Write-Host ""
Write-Host "Batch finished. Summary: $summaryPath"
