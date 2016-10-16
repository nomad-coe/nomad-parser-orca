package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object OrcaParser extends SimpleExternalParserGenerator(
  name = "OrcaParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("OrcaParser")) ::
      ("parserId" -> jn.JString("OrcaParser" + lab.OrcaVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.OrcaVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*
\s+\* O   R   C   A \*\s*
\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*
\s*
\s*--- An Ab Initio, DFT and Semiempirical electronic structure package ---\s*
""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/orca/parser/parser-orca/orca_parser.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-orca/orca_parser.py",
    "parser-orca/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/orca.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-orca" -> "parsers/orca/parser/parser-orca",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping(),
  metaInfoEnv = Some(lab.meta.KnownMetaInfoEnvs.orca)
)
