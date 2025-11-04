# Changelog

All notable changes to Orpheus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0-dev] - 2025-11-05

### Added
- ✨ **Integrated Report Generation**: New `integrated_report.tsv` file providing gene-level comprehensive information
  - Organized by gene with all information in one view
  - Includes BUSCO completeness status and gene ID
  - Shows multiple ORF completeness and length information
  - Displays homology annotation results (protein match, similarity, E-value)
  - Perfect for quick gene quality and function assessment
- ✨ **BUSCO Scoring System Enhancement**: Fragmented BUSCO genes now receive reasonable medium scores
  - Fragmented genes score 0.5 (previously 0.0)
  - Better reflects the partial functional evidence
  - More reasonable for downstream transcript selection
- ✨ **Intelligent File Discovery**: Auto-find intermediate files in working directory
  - No need to specify `-i` parameter when resuming from a specific step
  - Smart search for CD-HIT results, GFF3 files, and TransDecoder directories
  - Detailed documentation in pipeline flow guide

### Fixed
- 🐛 **BUSCO full_table.tsv Deep Directory Search**: Fixed issue where BUSCO results in nested subdirectories weren't found
  - Problem: Code only searched one level deep, but newer BUSCO creates deep structures like `busco_after/orpheus_after/run_lineage/full_table.tsv`
  - Solution: Use recursive glob pattern `**/full_table.tsv` to search all subdirectories
  - Impact: Resolved "BUSCO full_table.tsv not found" warnings in scoring step

### Changed
- 📝 Updated documentation with new integrated report features
- 🎯 Improved scoring system documentation with Fragmented BUSCO handling

## [0.1.0] - 2025-11-04

### Added
- ✨ **Comprehensive Scoring System**: Multi-dimensional transcript quality assessment
  - BUSCO completeness scoring (40% weight)
  - ORF completeness scoring (30% weight)
  - Homology evidence scoring (20% weight)
  - Sequence length scoring (10% weight)
- ✨ **Homology Search Support**: Integration with Diamond/BLASTP for protein database alignment
  - Improved ORF prediction accuracy
  - Support for SwissProt and custom databases
  - 100-10000x faster with Diamond
- ✨ **Complete ORF Filtering**: Identify complete ORFs with start and stop codons
  - Essential for PCR primer design
  - Configurable via `use_complete_only` parameter
- 📊 **Detailed TSV Reports**: 
  - `transcript_scores.tsv` - Per-transcript scoring details
  - Quality classification (high/medium/low)
- 🎯 **Flexible Pipeline Control**: Start from any step with `--start-from` parameter
- 🧵 **Multi-threading Support**: Thread count control via `-t` parameter
- 📝 **Help System**: Comprehensive help with `-h` and version info with `-v`

### Fixed
- 🐛 TransDecoder homology evidence integration
- 🐛 ORF completeness type parsing from GFF3

## [0.0.1] - 2025-11-01

### Added
- 🎉 Initial release
- CD-HIT redundancy removal
- TransDecoder ORF prediction
- Basic BUSCO integration
- YAML configuration system
- Command-line interface

---

## Legend

- ✨ New feature
- 🐛 Bug fix
- 📝 Documentation
- 🎯 Improvement
- 🔧 Configuration
- ⚠️ Deprecated
- 🗑️ Removed
- 🚀 Performance
- 🔒 Security

[Unreleased]: https://github.com/YanZezhang-debug/Orpheus/compare/v0.2.0-dev...HEAD
[0.2.0-dev]: https://github.com/YanZezhang-debug/Orpheus/compare/v0.1.0...v0.2.0-dev
[0.1.0]: https://github.com/YanZezhang-debug/Orpheus/compare/v0.0.1...v0.1.0
[0.0.1]: https://github.com/YanZezhang-debug/Orpheus/releases/tag/v0.0.1

