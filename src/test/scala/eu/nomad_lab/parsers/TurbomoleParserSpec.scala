package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object TurbomoleParserSpec extends Specification {
  "TurbomoleParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(TurbomoleParser, "parsers/turbomole/test/examples/NO_2_UHF.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(TurbomoleParser, "parsers/turbomole/test/examples/NO_2_UHF.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
