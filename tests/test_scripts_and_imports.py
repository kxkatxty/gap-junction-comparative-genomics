from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
TOOLS = ROOT / "tools"


class TestRepoLayoutSmoke:
    def test_scripts_exist_and_executable(self):
        runners = sorted(SCRIPTS.glob("run_*.sh"))
        assert len(runners) >= 10
        assert (SCRIPTS / "_env.sh").exists()
        for path in runners:
            mode = path.stat().st_mode
            assert mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH), f"not executable: {path.name}"

    def test_tools_python_modules_importable(self):
        modules = [
            "build_config_from_species",
            "build_quality_table",
            "classify_family_annotation",
            "collect_species_by_clade",
            "collect_unannotated_connexin_by_clade",
            "discover_innexins",
            "download_by_family",
            "download_gff3_from_species_list",
            "extract_exon_structures",
            "fetch_biocentral_predictions",
            "fetch_protein_features",
            "plot_exon_maps",
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{ROOT}:{TOOLS}:{env.get('PYTHONPATH', '')}"
        for name in modules:
            proc = subprocess.run(
                ["python3", "-c", f"import {name}"],
                cwd=ROOT,
                env=env,
                capture_output=True,
                text=True,
            )
            assert proc.returncode == 0, f"{name}: {proc.stderr}"

    def test_env_script_sets_pythonpath(self):
        proc = subprocess.run(
            ["bash", "-c", 'source scripts/_env.sh && python3 -c "import discover_innexins, pipeline.common"'],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr

    def test_pipeline_builders_importable(self):
        builders = [
            "pipeline.gap_junction_path_builder",
            "pipeline.inx_vs_cnx_builder",
            "pipeline.family_comparison_builder",
            "pipeline.innexin_similarity_builder",
            "pipeline.connexin_similarity_builder",
            "pipeline.innexin_clade_comparison_builder",
            "pipeline.connexin_clade_comparison_builder",
            "pipeline.new_species_gallery_builder",
            "pipeline.connexin_common",
            "pipeline.innexin_reference_panel_builder",
            "pipeline.merge_curator_probe_into_master",
            "pipeline.batch_innexin_curator_probe",
            "pipeline.gff_synteny_neighbor_summary",
            "pipeline.innexin_subfamily_hypothesis",
            "pipeline.phylogenetic_story_builder",
            "pipeline.synteny_case_study_1_dmel_cluster",
            "pipeline.synteny_phylo_guided_inx_cluster",
        ]
        for mod in builders:
            proc = subprocess.run(
                ["python3", "-c", f"import {mod}"],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            assert proc.returncode == 0, f"{mod}: {proc.stderr}"
