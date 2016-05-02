package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object OrcaParserSpec extends Specification {
  "OrcaParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OrcaParser, "parsers/orca/test/examples/", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OrcaParser, "parsers/orca/test/examples/", "json") must_== ParseResult.ParseSuccess
    }
  }
}
