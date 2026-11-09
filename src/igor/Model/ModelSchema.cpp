// ModelSchema.cpp ---

#include <igor/Model/ModelSchema.h>

#include <nlohmann/json-schema.hpp>
#include <nlohmann/json.hpp>

#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace igor::model::schema {

namespace detail {
/// Generated from the schema files by cmake/igorEmbedSchemas.cmake.
const std::map<std::string, std::string_view> &embedded_schemas();
}  // namespace detail

namespace {

constexpr const char *kCombinationSummary = "at least one subschema has failed";

/// Every violation, at its JSON pointer: one message lists them all.
class CollectingErrorHandler : public nlohmann::json_schema::basic_error_handler
{
public:
    std::vector<std::string> errors;

    void error(const nlohmann::json::json_pointer &pointer, const nlohmann::json &instance,
               const std::string &message) override
    {
        basic_error_handler::error(pointer, instance, message);
        // An allOf that fails reports a summary and then its case, each nested allOf adding its
        // prefix: keep the case, without the prefixes. The tagged unions are allOf's of if/then.
        std::string text = message;
        while (text.rfind("[combination: ", 0) == 0)
            text = text.substr(text.find("] ") + 2);
        if (text.rfind(kCombinationSummary, 0) == 0)
            return;
        // additionalProperties: false, in the words of a model file's reader.
        const std::string additional = "validation failed for additional property '";
        const std::string refused = "': instance invalid as per false-schema";
        if (text.rfind(additional, 0) == 0 and text.size() > additional.size() + refused.size()
            and text.compare(text.size() - refused.size(), refused.size(), refused) == 0)
            text = "unknown key \"" + text.substr(additional.size(), text.size() - additional.size() - refused.size())
                   + "\"";
        std::string entry = (pointer.empty() ? std::string("/") : pointer.to_string()) + ": " + text;
        // A value says what was wrong with it; an object or a list would only make the line long.
        if (instance.is_primitive())
            entry += " (" + instance.dump() + ")";
        errors.push_back(std::move(entry));
    }
};

/// Serves a $ref to an embedded schema, by the file name that ends its URI.
void load_embedded(const nlohmann::json_uri &uri, nlohmann::json &schema)
{
    const std::string &path = uri.path();
    const std::string file = path.substr(path.find_last_of('/') + 1);
    schema = nlohmann::json::parse(schema_text(file));
}

/// The validator of one schema file. Built once.
const nlohmann::json_schema::json_validator &validator_for(const std::string &file)
{
    static std::mutex mutex;
    static std::map<std::string, std::unique_ptr<nlohmann::json_schema::json_validator>> validators;

    const std::lock_guard<std::mutex> lock(mutex);
    auto &slot = validators[file];
    if (not slot) {
        nlohmann::json root;
        root["$ref"] = std::string(kBaseUri) + file;
        slot = std::make_unique<nlohmann::json_schema::json_validator>(root, load_embedded);
    }
    return *slot;
}

}  // namespace

std::string document_schema_uri()
{
    return std::string(kBaseUri) + "model.schema.json";
}

std::vector<std::string> schema_files()
{
    std::vector<std::string> out;
    for (const auto &entry : detail::embedded_schemas())
        out.push_back(entry.first);
    return out;
}

std::string_view schema_text(const std::string &file)
{
    const auto &files = detail::embedded_schemas();
    const auto found = files.find(file);
    if (found == files.end())
        throw std::out_of_range("model schema: no schema file \"" + file + "\"");
    return found->second;
}

void validate(const nlohmann::json &instance, const std::string &file, const std::string &what,
              std::initializer_list<const char *> required)
{
    CollectingErrorHandler handler;
    validator_for(file).validate(instance, handler);
    // What a reader needs on top of the schema, which leaves every section optional.
    for (const char *name : required) {
        if (instance.is_object() and not instance.contains(name))
            handler.errors.push_back(std::string("/: required property '") + name + "' not found in object");
    }
    if (handler.errors.empty())
        return;

    // A document can break the schema in thousands of places; the first few say what is wrong.
    constexpr std::size_t kShown = 10;
    std::string message = "model json: " + what + " does not follow " + file + ":";
    for (std::size_t i = 0; i != handler.errors.size() and i != kShown; ++i)
        message += "\n  " + handler.errors[i];
    if (handler.errors.size() > kShown)
        message += "\n  and " + std::to_string(handler.errors.size() - kShown) + " more";
    throw std::runtime_error(message);
}

}  // namespace igor::model::schema

//
// ModelSchema.cpp ends here
