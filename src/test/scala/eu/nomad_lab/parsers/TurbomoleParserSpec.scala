package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object TurbomoleParserSpec extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_bandstructure/B3LYP/Si2.castep_b_sp_v_1", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_bandstructure/B3LYP/Si2.castep_b_sp_v_1", "json") must_== ParseResult.ParseSuccess
    }
  }
}
